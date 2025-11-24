/* ----------------------------------------------------------
   0. (Optional) set schema
   ---------------------------------------------------------- */
-- USE your_schema_name;


/* ----------------------------------------------------------
   1. Set variance thresholds (balanced & reliable)
   ---------------------------------------------------------- */

SET @thr_gap_min   = 0.03;   -- stable gap
SET @thr_gap_alpha = 0.03;
SET @thr_gap_beta  = 0.03;

SET @thr_dipmom    = 0.15;   -- conformer-level noise allowed

-- How much we trust large gap (onset_norm) vs weak absorption
SET @w_onset       = 0.7;   -- weight for onset_norm
SET @w_absorption  = 0.3;   -- weight for (1 - absorption_strength_norm)

-- Min PCE for "useful" OPVs (optional)
SET @pce_min       = 3.0;   -- tweak later

/* ----------------------------------------------------------
   2. Per-molecule stats: averages & variances by mol_graph_id
      for data_calcqcset1
   ---------------------------------------------------------- */

DROP VIEW IF EXISTS v_calcqcset1_stats;

CREATE VIEW v_calcqcset1_stats AS
SELECT
    mol_graph_id,
    COUNT(*)                AS n_experiments,

    AVG(e_gap_min)          AS avg_e_gap_min,
    VAR_POP(e_gap_min)      AS var_e_gap_min,

    AVG(e_gap_alpha)        AS avg_e_gap_alpha,
    VAR_POP(e_gap_alpha)    AS var_e_gap_alpha,

    AVG(e_gap_beta)         AS avg_e_gap_beta,
    VAR_POP(e_gap_beta)     AS var_e_gap_beta,

    AVG(dipmom_total)       AS avg_dipmom_total,
    VAR_POP(dipmom_total)   AS var_dipmom_total
FROM data_calcqcset1
GROUP BY mol_graph_id;


/* ----------------------------------------------------------
   3. Clean table for data_calcqcset1:
      - One row per mol_graph_id
      - Keep only low-variance molecules
      - Store averages + variances + count
   ---------------------------------------------------------- */

DROP TABLE IF EXISTS data_calcqcset1_clean;

CREATE TABLE data_calcqcset1_clean AS
SELECT
    MIN(id)                      AS rep_calc_id,           -- representative calc row
    mol_graph_id,
    MIN(mol_geom_id)             AS rep_mol_geom_id,       -- representative geometry
    COUNT(*)                     AS n_experiments,

    AVG(e_gap_min)               AS e_gap_min_avg,
    VAR_POP(e_gap_min)           AS e_gap_min_var,

    AVG(e_gap_alpha)             AS e_gap_alpha_avg,
    VAR_POP(e_gap_alpha)         AS e_gap_alpha_var,

    AVG(e_gap_beta)              AS e_gap_beta_avg,
    VAR_POP(e_gap_beta)          AS e_gap_beta_var,

    AVG(dipmom_total)            AS dipmom_total_avg,
    VAR_POP(dipmom_total)        AS dipmom_total_var
FROM data_calcqcset1
GROUP BY mol_graph_id
HAVING
    -- Treat NULL variance (e.g., single experiment or all NULLs) as acceptable
    (VAR_POP(e_gap_min)    <= @thr_gap_min   OR VAR_POP(e_gap_min)    IS NULL)
AND (VAR_POP(e_gap_alpha)  <= @thr_gap_alpha OR VAR_POP(e_gap_alpha)  IS NULL)
AND (VAR_POP(e_gap_beta)   <= @thr_gap_beta  OR VAR_POP(e_gap_beta)   IS NULL)
AND (VAR_POP(dipmom_total) <= @thr_dipmom    OR VAR_POP(dipmom_total) IS NULL);


/* ----------------------------------------------------------
   4. Optional: which mol_graph_ids were dropped (calcqc)?
   ---------------------------------------------------------- */

DROP VIEW IF EXISTS v_calcqcset1_dropped;

CREATE VIEW v_calcqcset1_dropped AS
SELECT
    s.*
FROM v_calcqcset1_stats s
LEFT JOIN data_calcqcset1_clean c
  ON s.mol_graph_id = c.mol_graph_id
WHERE c.mol_graph_id IS NULL;


/* ----------------------------------------------------------
   5. Per-molecule stats for data_calibqcset1
      (grouped by mol_graph_id + calib_type)
   ---------------------------------------------------------- */

DROP VIEW IF EXISTS v_calibqcset1_stats;

CREATE VIEW v_calibqcset1_stats AS
SELECT
    mol_graph_id,
    calib_type,
    COUNT(*)              AS n_experiments,

    AVG(e_homo_alpha)     AS avg_e_homo_alpha,
    VAR_POP(e_homo_alpha) AS var_e_homo_alpha,

    AVG(e_lumo_alpha)     AS avg_e_lumo_alpha,
    VAR_POP(e_lumo_alpha) AS var_e_lumo_alpha,

    AVG(e_gap_alpha)      AS avg_e_gap_alpha,
    VAR_POP(e_gap_alpha)  AS var_e_gap_alpha,

    AVG(e_homo_beta)      AS avg_e_homo_beta,
    VAR_POP(e_homo_beta)  AS var_e_homo_beta,

    AVG(e_lumo_beta)      AS avg_e_lumo_beta,
    VAR_POP(e_lumo_beta)  AS var_e_lumo_beta,

    AVG(e_gap_beta)       AS avg_e_gap_beta,
    VAR_POP(e_gap_beta)   AS var_e_gap_beta,

    AVG(e_gap_min)        AS avg_e_gap_min,
    VAR_POP(e_gap_min)    AS var_e_gap_min
FROM data_calibqcset1
GROUP BY mol_graph_id, calib_type;


/* ----------------------------------------------------------
   6. Clean table for data_calibqcset1:
      - One row per (mol_graph_id, calib_type)
      - Keep only low-variance combos
      - Store averages + variances + count
   ---------------------------------------------------------- */

DROP TABLE IF EXISTS data_calibqcset1_clean;

CREATE TABLE data_calibqcset1_clean AS
SELECT
    MIN(id)                 AS rep_calib_id,         -- representative row id
    mol_graph_id,
    MIN(mol_geom_id)        AS rep_mol_geom_id,      -- representative geometry
    calib_type,
    COUNT(*)                AS n_experiments,

    AVG(e_homo_alpha)       AS e_homo_alpha_avg,
    VAR_POP(e_homo_alpha)   AS e_homo_alpha_var,

    AVG(e_lumo_alpha)       AS e_lumo_alpha_avg,
    VAR_POP(e_lumo_alpha)   AS e_lumo_alpha_var,

    AVG(e_gap_alpha)        AS e_gap_alpha_avg,
    VAR_POP(e_gap_alpha)    AS e_gap_alpha_var,

    AVG(e_homo_beta)        AS e_homo_beta_avg,
    VAR_POP(e_homo_beta)    AS e_homo_beta_var,

    AVG(e_lumo_beta)        AS e_lumo_beta_avg,
    VAR_POP(e_lumo_beta)    AS e_lumo_beta_var,

    AVG(e_gap_beta)         AS e_gap_beta_avg,
    VAR_POP(e_gap_beta)     AS e_gap_beta_var,

    AVG(e_gap_min)          AS e_gap_min_avg,
    VAR_POP(e_gap_min)      AS e_gap_min_var
FROM data_calibqcset1
GROUP BY mol_graph_id, calib_type
HAVING
    -- Use the same gap-variance thresholds as before.
    -- Treat NULL variance (e.g. only one experiment) as acceptable.
    (VAR_POP(e_gap_min)   <= @thr_gap_min   OR VAR_POP(e_gap_min)   IS NULL)
AND (VAR_POP(e_gap_alpha) <= @thr_gap_alpha OR VAR_POP(e_gap_alpha) IS NULL)
AND (VAR_POP(e_gap_beta)  <= @thr_gap_beta  OR VAR_POP(e_gap_beta)  IS NULL);


/* ----------------------------------------------------------
   7. View of dropped (high-variance) mol_graph_id + calib_type
   ---------------------------------------------------------- */

DROP VIEW IF EXISTS v_calibqcset1_dropped;

CREATE VIEW v_calibqcset1_dropped AS
SELECT
    s.*
FROM v_calibqcset1_stats s
LEFT JOIN data_calibqcset1_clean c
  ON  s.mol_graph_id = c.mol_graph_id
  AND s.calib_type   = c.calib_type
WHERE c.mol_graph_id IS NULL;


/* ----------------------------------------------------------
   8. Aggregate Scharber outputs per mol_graph_id
      (in case there are multiple rows per molecule)
   ---------------------------------------------------------- */

CREATE OR REPLACE VIEW v_scharber_agg AS
SELECT
    mol_graph_id,
    AVG(jsc) AS jsc_avg,
    AVG(pce) AS pce_avg
FROM data_scharber
GROUP BY mol_graph_id;


/* ----------------------------------------------------------
   9. Absorption-related features + transparency proxy
   ---------------------------------------------------------- */

CREATE OR REPLACE VIEW v_absorption_features AS
WITH base AS (
  SELECT
    g.id           AS mol_graph_id,
    g.SMILES_str,
    g.n_heavyatoms,
    g.mass,

    -- Cleaned gaps
    c.e_gap_min_avg   AS e_gap_min,
    c.e_gap_alpha_avg AS e_gap_alpha,
    c.e_gap_beta_avg  AS e_gap_beta,

    -- Scharber metrics
    s.jsc_avg AS jsc,
    s.pce_avg AS pce,

    -- Spin asymmetry
    ABS(c.e_gap_alpha_avg - c.e_gap_beta_avg) AS gap_asymmetry
  FROM data_molgraph g
  JOIN data_calcqcset1_clean c ON g.id = c.mol_graph_id
  JOIN v_scharber_agg        s ON g.id = s.mol_graph_id
  WHERE
    -- require Scharber
    s.pce_avg BETWEEN 0 AND 40
    AND s.jsc_avg BETWEEN 0 AND 30
    -- no badly broken spin symmetry
    AND ABS(c.e_gap_alpha_avg - c.e_gap_beta_avg) < 0.2
    -- optional: only somewhat functional devices
    AND s.pce_avg >= @pce_min
)
SELECT
  b.*,

  -- Normalized Jsc across this OPV-like set
  (b.jsc - MIN(b.jsc) OVER()) /
  NULLIF(MAX(b.jsc) OVER() - MIN(b.jsc) OVER(), 0) AS absorption_strength_norm,

  -- Normalized onset (0..1)
  (b.e_gap_min - MIN(b.e_gap_min) OVER()) /
  NULLIF(MAX(b.e_gap_min) OVER() - MIN(b.e_gap_min) OVER(), 0) AS onset_norm,

  -- Transparency proxy (raw, before normalization):
  -- high when onset high AND absorption_strength low
  (
      @w_onset      * (
          (b.e_gap_min - MIN(b.e_gap_min) OVER()) /
          NULLIF(MAX(b.e_gap_min) OVER() - MIN(b.e_gap_min) OVER(), 0)
      )
    + @w_absorption * (
          1 - (
              (b.jsc - MIN(b.jsc) OVER()) /
              NULLIF(MAX(b.jsc) OVER() - MIN(b.jsc) OVER(), 0)
          )
      )
  ) AS transparency_raw

FROM base b;

/* ----------------------------------------------------------
   10. Final dataset view for GNN
       - Regression: transparency_pct (0..100)
       - Plus raw/z-scored targets for experiments
   ---------------------------------------------------------- */

CREATE OR REPLACE VIEW v_gnn_dataset AS
WITH t AS (
  SELECT
    v.*,
    -- Normalize transparency_raw into [0,1]
    (transparency_raw - MIN(transparency_raw) OVER()) /
    NULLIF(MAX(transparency_raw) OVER() - MIN(transparency_raw) OVER(), 0) AS transparency_norm
  FROM v_absorption_features v
)
SELECT
    mol_graph_id,
    SMILES_str,
    n_heavyatoms,
    mass,
    e_gap_min,
    e_gap_alpha,
    e_gap_beta,
    jsc,
    pce,
    gap_asymmetry,
    absorption_strength_norm,
    onset_norm,
    transparency_raw,
    transparency_norm,
    transparency_norm * 100.0 AS transparency_pct,  -- << YOUR main regression target

    -- Optional binary label: "high transparency"
    CASE
        WHEN transparency_norm >= 0.75 THEN 1   -- top 25% transparent
        ELSE 0
    END AS y_transparent_high

FROM t;
-- End of pipeline.sql