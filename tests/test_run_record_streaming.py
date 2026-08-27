from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUN_RECORD = REPO / "tools" / "run_record.py"


def test_command_output_reaches_the_log_before_the_command_finishes(
    tmp_path: Path,
) -> None:
    """A multi-hour build is only monitorable if its log streams.

    ``read(n)`` on a blocking pipe waits for the full n bytes, so output used to
    surface in 64 KiB steps and the pending tail was lost outright if the
    process was killed -- which made the log actively misleading about progress.
    """
    record = tmp_path / "probe"
    subprocess.run(
        [sys.executable, str(RUN_RECORD), "init", str(record),
         "--kind", "experiment", "--title", "probe"],
        cwd=REPO, check=True, capture_output=True,
    )
    emitter = tmp_path / "emit.py"
    emitter.write_text(
        "import time\n"
        "print('EARLY', flush=True)\n"
        "time.sleep(8)\n"
        "print('LATE', flush=True)\n"
    )

    process = subprocess.Popen(
        [sys.executable, str(RUN_RECORD), "exec", str(record),
         "--label", "probe", "--", sys.executable, str(emitter)],
        cwd=REPO, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    log = record / "logs" / "probe.log"
    try:
        deadline = time.monotonic() + 5
        seen = ""
        while time.monotonic() < deadline:
            if log.exists():
                seen = log.read_text()
                if "EARLY" in seen:
                    break
            time.sleep(0.2)
        assert "EARLY" in seen, (
            "early output never reached the log while the command was still "
            f"running (log held {seen!r})"
        )
        assert "LATE" not in seen, "sanity: the command should still be running"
    finally:
        process.wait(timeout=30)

    assert log.read_text().split() == ["EARLY", "LATE"]
