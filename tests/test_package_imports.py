"""The package's module graph stays acyclic, and every module imports alone.

A cycle only fails in one import order. Importing ``cli`` pulls in nearly
everything in a fixed sequence, so a reintroduced ``a -> b -> a`` can sit
undetected behind whichever module the suite happens to load first. Each module
therefore gets its own interpreter here, entered directly.

This guards the arrangement the ``sources``/``build_store``/``build_lock`` split
exists to create: ``sources`` holds the declarative model, ``build_store`` the
engine, and ``build_lock`` the record -- so the engine and the record can both
depend on the model without depending on each other.
"""

from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import pytest

import gks_refgetstore

PACKAGE = Path(gks_refgetstore.__file__).parent
MODULES = sorted(p.stem for p in PACKAGE.glob("*.py") if p.stem != "__init__")


def test_the_module_list_is_not_silently_empty() -> None:
    """A bad glob would make every parametrized case vanish rather than fail."""
    assert len(MODULES) >= 12


@pytest.mark.parametrize("module", MODULES)
def test_each_module_imports_in_isolation(module: str) -> None:
    """A fresh interpreter per module, so import order cannot mask a cycle."""
    result = subprocess.run(
        [sys.executable, "-c", f"import gks_refgetstore.{module}"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, (
        f"gks_refgetstore.{module} failed to import on its own:\n{result.stderr}"
    )


def _intra_package_imports(path: Path) -> set[str]:
    """Module-level first-party dependencies, read statically.

    Static rather than runtime because a cycle broken by deferring an import
    into a function body is still a cycle in the design; the point of the graph
    is that no module *needs* deferring. Imports nested inside a function are
    skipped for exactly that reason -- the remaining ones are deliberate
    (``gtars``, kept out of module scope so ``--help`` stays cheap).
    """
    tree = ast.parse(path.read_text())
    deps: set[str] = set()
    for node in tree.body:  # module level only
        if isinstance(node, ast.ImportFrom) and node.level == 1:
            # ``from . import a, b`` names modules; ``from .a import x`` does not
            deps |= ({alias.name for alias in node.names} if node.module is None
                     else {node.module})
    return deps & set(MODULES)


def test_the_module_graph_is_acyclic() -> None:
    graph = {m: _intra_package_imports(PACKAGE / f"{m}.py") for m in MODULES}

    def find_cycle(node: str, stack: list[str]) -> list[str] | None:
        if node in stack:
            return stack[stack.index(node):] + [node]
        for dep in sorted(graph[node]):
            if (cycle := find_cycle(dep, stack + [node])) is not None:
                return cycle
        return None

    cycles = [c for m in MODULES if (c := find_cycle(m, [])) is not None]
    assert not cycles, "import cycle: " + " -> ".join(cycles[0])


def test_importing_the_package_does_not_pull_in_the_engine() -> None:
    """``import gks_refgetstore`` must stay cheap.

    ``__init__`` re-exporting submodules would drag in gtars and the whole
    ingest path, which is what the CLI's lazy dispatch exists to avoid.
    """
    result = subprocess.run(
        [sys.executable, "-c",
         "import sys, gks_refgetstore; "
         "print([m for m in sys.modules if m.startswith('gks_refgetstore.') "
         "or m == 'gtars'])"],
        capture_output=True, text=True, check=True,
    )
    assert result.stdout.strip() == "[]", result.stdout
