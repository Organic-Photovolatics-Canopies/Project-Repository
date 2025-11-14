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
   9. Absorption-related features using CLEANED calc table
      + cleaning filters
   ---------------------------------------------------------- */

CREATE OR REPLACE VIEW v_absorption_features AS
SELECT
  g.id AS mol_graph_id,
  g.SMILES_str,
  g.n_heavyatoms,
  g.mass,

  -- Use averaged, cleaned gaps from data_calcqcset1_clean
  c.e_gap_min_avg   AS e_gap_min,
  c.e_gap_alpha_avg AS e_gap_alpha,
  c.e_gap_beta_avg  AS e_gap_beta,

  -- Aggregated Scharber metrics
  s.jsc_avg AS jsc,
  s.pce_avg AS pce,

  -- Derived features
  ABS(c.e_gap_alpha_avg - c.e_gap_beta_avg) AS gap_asymmetry,

  -- Normalized Jsc across all molecules with Jsc
  (s.jsc_avg - MIN(s.jsc_avg) OVER()) /
  NULLIF(MAX(s.jsc_avg) OVER() - MIN(s.jsc_avg) OVER(), 0) AS absorption_strength_norm,

  -- Normalized onset (gap_min) across all molecules with gaps
  (c.e_gap_min_avg - MIN(c.e_gap_min_avg) OVER()) /
  NULLIF(MAX(c.e_gap_min_avg) OVER() - MIN(c.e_gap_min_avg) OVER(), 0) AS onset_norm

FROM data_molgraph g
JOIN data_calcqcset1_clean c
       ON g.id = c.mol_graph_id
JOIN v_scharber_agg s
       ON g.id = s.mol_graph_id
WHERE
  -- require valid Scharber entry
  s.mol_graph_id IS NOT NULL
  -- physically reasonable Scharber outputs
  AND s.pce_avg BETWEEN 0 AND 40
  AND s.jsc_avg BETWEEN 0 AND 30
  -- no badly broken spin symmetry
  AND ABS(c.e_gap_alpha_avg - c.e_gap_beta_avg) < 0.2;


/* ----------------------------------------------------------
   10. Final GNN dataset view with binary transparency label
   ---------------------------------------------------------- */

CREATE OR REPLACE VIEW v_gnn_dataset AS
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
    CASE
        WHEN onset_norm >= 0.5 THEN 1  -- ≥ 50% transparent
        ELSE 0
    END AS y_transparent
FROM v_absorption_features;
-- End of pipeline.sql