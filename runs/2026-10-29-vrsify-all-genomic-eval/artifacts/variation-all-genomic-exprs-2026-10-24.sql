-------------------------------------------------------------------------------
-- variation_all_genomic_exprs — every genomic expression of every ClinVar variation
--
-- A one-off companion to variation_identity for measuring VRS translation at scale,
-- not a pipeline step. variation_identity keeps ONE expression per variation (the
-- highest-precedence candidate); this keeps EVERY genomic candidate that
-- variation_identity_build's Step 7 chooses among, without its per-accession and
-- per-variation dedup, plus every alternate HGVS spelling ClinVar gives on the same
-- accession:
--
--   precedence 1  spdi    canonical SPDI
--   precedence 2  hgvs    genomic, top-level HGVS (every assembly)
--   precedence 3  gnomad  gnomAD-style id from the location (every assembly)
--   precedence 4  hgvs    HGVS derived from the location, where no gnomAD id exists
--   precedence 5  hgvs    genomic HGVS that is not top-level (NG_, alternate loci, patches)
--
-- An expression string that two ranks produce is kept once, at the higher rank.
-- So one allele appears once per format and per assembly, and the VRS ids of the
-- expressions on the same sequence can be compared with each other.
--
-- Reads only tables variation_identity_build already persisted for the release
-- (variation_identity, variation_hgvs, variation_loc); it does not modify them.
-- Rows are limited to the VRS classes vrsify translates. Columns follow
-- variation_identity's, so vrsify reads the export unchanged, plus:
--   hgvs_type           variation_hgvs.type of an hgvs row
--   clinvar_expression  the expression as ClinVar gives it, before the del<N> cleanup
--   selected            this row is the one variation_identity selected
--
-- Writes `<release schema>.variation_all_genomic_exprs` and exports it as sharded
-- newline-delimited JSON to export_uri, which must contain one '*'. The export never
-- overwrites: it fails if objects already match export_uri.
--
--   CALL `clinvar_ingest.variation_all_genomic_exprs`(DATE '2026-10-24',
--     'gs://clinvar-gkm/2026-10-24/dev/scratch/all-genomic/vi-all-genomic-exprs-*.jsonl.gz');
-------------------------------------------------------------------------------
DECLARE on_date DATE DEFAULT DATE '2026-10-24';
DECLARE export_uri STRING DEFAULT 'gs://clinvar-gkm/2026-10-24/dev/scratch/all-genomic/vi-all-genomic-exprs-*.jsonl.gz';
BEGIN
  DECLARE query_all_genomic_exprs STRING;
  DECLARE query_export STRING;

  IF NOT REGEXP_CONTAINS(export_uri, r'^gs://[^*]+\*[^*]*\.jsonl\.gz$') THEN
    RAISE USING MESSAGE = FORMAT('export_uri must be a gs:// URI with one * ending in .jsonl.gz: %s', export_uri);
  END IF;

  FOR rec IN (select s.schema_name FROM clinvar_ingest.schema_on(on_date) as s)
  DO

    -------------------------------------------------------------------------
    -- Collect the genomic candidates (Step 7's ranks 1-5) without dedup
    -------------------------------------------------------------------------
    SET query_all_genomic_exprs = REPLACE("""
      CREATE OR REPLACE TABLE `{S}.variation_all_genomic_exprs` AS
      WITH vi AS (
        select
          variation_id,
          vrs_class,
          absolute_copies,
          range_copies,
          variation_type,
          subclass_type,
          assembly_version as vi_assembly_version,
          fmt as vi_fmt,
          source as vi_source,
          precedence as vi_precedence
        from `{S}.variation_identity`
        where
          vrs_class in ('Allele', 'CopyNumberChange', 'CopyNumberCount')
      ),
      h AS (
        -- every nucleotide spelling on each genomic accession: variation_hgvs keeps one
        -- representative per accession in hgvs_source and all of them in expr
        select
          vh.variation_id,
          vh.assembly_version,
          vh.accession,
          vh.type as hgvs_type,
          REGEXP_REPLACE(e.nucleotide, r"del[0-9]+", "del") as source,
          MIN(e.nucleotide) as clinvar_expression,
          LOGICAL_AND(MOD(LENGTH(REGEXP_REPLACE(e.nucleotide, r"[^\\(\\)]", "")), 2) = 0) AS has_balanced_parens,
          IF(vh.type = 'genomic, top-level', 2, 5) as precedence
        from `{S}.variation_hgvs` vh
        cross join unnest(vh.expr) as e
        where
          vh.type in ('genomic, top-level', 'genomic')
          and
          e.nucleotide is not null
        group by
          vh.variation_id,
          vh.assembly_version,
          vh.accession,
          vh.type,
          source
      ),
      h_issues AS (
        -- variation_identity_build Step 3's issues, applied to every spelling rather
        -- than only the representative
        SELECT
          h.*,
          CASE
          WHEN (NOT REGEXP_CONTAINS(h.source, r'^(NC|NT|NW|NG|NM|NR|XM|XR)_')) THEN
            'sequence for accession not supported by vrs-python release'
          WHEN REGEXP_CONTAINS(h.source, r'\\[[\\(\\)\\-0-9]+\\]') THEN
            'repeat expressions are not supported.'
          WHEN NOT h.has_balanced_parens THEN
            'expression contains unbalaned paretheses.'
          WHEN REGEXP_CONTAINS(h.source, r'[0-9]+[\\+\\-][0-9]+') THEN
            'intronic positions are not resolvable in sequence.'
          WHEN REGEXP_CONTAINS(h.source, r'^NP') THEN
            'protein expressions not supported.'
          WHEN NOT (
            --snv
            REGEXP_CONTAINS(h.source, r'^(NC|NT|NW|NG|NM|NR|XM|XR)_[0-9]+\\.[0-9]+\\:[gmcnr]\\.[0-9]+[ACTGN]\\>[ACTGN]+$') OR
            -- same as ref
            REGEXP_CONTAINS(h.source, r'^(NC|NT|NW|NG|NM|NR|XM|XR)_[0-9]+\\.[0-9]+\\:[gmcnr]\\.[0-9]+[ACTGN]?\\=$') OR
            -- single residue dup or del or delins?
            REGEXP_CONTAINS(h.source, r'^(NC|NT|NW|NG|NM|NR|XM|XR)_[0-9]+\\.[0-9]+\\:[gmcnr]\\.[0-9]+(dup|del|delins)[ACTGN]*$') OR
            -- precise range dup or del or delins or ins
            REGEXP_CONTAINS(h.source, r'^(NC|NT|NW|NG|NM|NR|XM|XR)_[0-9]+\\.[0-9]+\\:[gmcnr]\\.[0-9]+_[0-9]+(dup|del|delins|ins)[ACTGN]*$') OR
            -- inner/outer range dup or del
            REGEXP_CONTAINS(h.source, r'^(NC|NT|NW|NG|NM|NR|XM|XR)_[0-9]+\\.[0-9]+\\:[gmcnr]\\.\\([0-9\\?]+_[0-9\\?]+\\)_\\([0-9\\?]+_[0-9\\?]+\\)(dup|del)[ACTGN]*$')
          ) THEN
            'unsupported hgvs expression.'
          END as issue
        from h
      ),
      candidates AS (
        select
          variation_id,
          vi_assembly_version as assembly_version,
          SPLIT(vi_source, ':')[OFFSET(0)] as accession,
          'spdi' as fmt,
          vi_source as source,
          CAST(null AS STRING) as issue,
          -- #1 spdi; always selected when present, so variation_identity holds all of them
          1 as precedence,
          CAST(null AS STRING) as hgvs_type,
          vi_source as clinvar_expression
        from vi
        where vi_fmt = 'spdi'
        UNION ALL
        select
          variation_id,
          assembly_version,
          accession,
          'hgvs' as fmt,
          source,
          issue,
          -- #2 hgvs genomic, top-level; #5 hgvs genomic (not top-level)
          precedence,
          hgvs_type,
          clinvar_expression
        from h_issues
        UNION ALL
        select DISTINCT
          variation_id,
          assembly_version,
          accession,
          'gnomad' as fmt,
          gnomad_source as source,
          CAST(null AS STRING) as issue,
          -- #3 gnomad location-based (genomic 'top-level')
          3 as precedence,
          CAST(null AS STRING) as hgvs_type,
          gnomad_source as clinvar_expression
        from `{S}.variation_loc`
        where
          gnomad_source is not null
        UNION ALL
        select DISTINCT
          variation_id,
          assembly_version,
          accession,
          'hgvs' as fmt,
          loc_hgvs_source as source,
          loc_hgvs_issue as issue,
          -- #4 derived hgvs for non-precise location regions (genomic 'top-level')
          4 as precedence,
          CAST(null AS STRING) as hgvs_type,
          loc_hgvs_source as clinvar_expression
        from `{S}.variation_loc`
        where
          loc_hgvs_source is not null
          and
          gnomad_source is null
      ),
      distinct_candidates AS (
        -- a location-derived hgvs (#4) can be the same string as ClinVar's own (#2);
        -- keep each expression once, at its highest precedence
        select
          *
        from candidates
        QUALIFY ROW_NUMBER() OVER (
          PARTITION BY variation_id, assembly_version, fmt, source
          ORDER BY precedence
        ) = 1
      )
      select
        c.variation_id,
        c.assembly_version,
        c.accession,
        vi.vrs_class,
        vi.absolute_copies,
        vi.range_copies,
        c.fmt,
        c.source,
        IF(
          vi.vrs_class = 'CopyNumberChange',
          CASE
            WHEN vi.variation_type IN ('Deletion', 'copy number loss') THEN
              "loss"
            WHEN vi.variation_type IN ('Duplication', 'copy number gain') THEN
              "gain"
            ELSE
              NULL
            END,
          NULL
        ) as copy_change_type,
        c.issue,
        c.precedence,
        vi.variation_type,
        vi.subclass_type,
        c.hgvs_type,
        c.clinvar_expression,
        IFNULL(
          c.precedence = vi.vi_precedence
          and c.fmt = vi.vi_fmt
          and c.source = vi.vi_source
          and IFNULL(c.assembly_version, 0) = IFNULL(vi.vi_assembly_version, 0),
          FALSE
        ) as selected
      from distinct_candidates c
      join vi
      on
        vi.variation_id = c.variation_id
    """, '{S}', rec.schema_name);

    EXECUTE IMMEDIATE query_all_genomic_exprs;

    -------------------------------------------------------------------------
    -- Export with integers as JSON strings, as `bq extract` writes the
    -- variation_identity export; vrsify looks assemblies up by "38", "37", "36"
    -------------------------------------------------------------------------
    SET query_export = REPLACE(REPLACE("""
      EXPORT DATA OPTIONS (
        uri = '{U}',
        format = 'JSON',
        compression = 'GZIP',
        overwrite = false
      ) AS
      select
        variation_id,
        CAST(assembly_version AS STRING) as assembly_version,
        accession,
        vrs_class,
        CAST(absolute_copies AS STRING) as absolute_copies,
        ARRAY(select CAST(rc AS STRING) from unnest(range_copies) as rc) as range_copies,
        fmt,
        source,
        copy_change_type,
        issue,
        CAST(precedence AS STRING) as precedence,
        variation_type,
        subclass_type,
        hgvs_type,
        clinvar_expression,
        selected
      from `{S}.variation_all_genomic_exprs`
    """, '{S}', rec.schema_name), '{U}', export_uri);

    EXECUTE IMMEDIATE query_export;

  END FOR;
END;

