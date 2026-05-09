# Follow-Up Priorities

## What this report means

This report does **not** identify damaged genes, DNA lesions, cancer, or clinical injury. It prioritizes RNA-expression signals in DNA damage-response pathways for **follow-up** under this repository’s exploratory scoring rules.

## Highest-priority pathway signals

Source: `outputs/tables/main_candidate_signals.csv` (main candidate gate applied by this pipeline).

### OSD-570_PBMC_snRNA | Other | DNA_repair

- **Comparison:** (R+45)_vs_(R+1)
- **pathway_rank_metric:** 1.619363931094595
- **average_signed_effect:** 0.4135154914347663
- **overlapping_gene_count:** 3
- **Why it deserves attention:** Exploratory pathway ranking with overlap ≥3 and main-candidate filters passed; interpretation label: moderate positive DNA damage-response signature.
- **Caution:** limited pathway-level signal; small crew sizes possible; not clinical severity.

### OSD-570_PBMC_snRNA | CD16+ Monocyte | oxidative_stress

- **Comparison:** (R+45)_vs_(R+1)
- **pathway_rank_metric:** 1.135604041940302
- **average_signed_effect:** -0.3466468461270288
- **overlapping_gene_count:** 4
- **Why it deserves attention:** Exploratory pathway ranking with overlap ≥3 and main-candidate filters passed; interpretation label: pathway genes changed but mostly downregulated.
- **Caution:** limited pathway-level signal; small crew sizes possible; not clinical severity.

### OSD-570_PBMC_snRNA | CD16+ Monocyte | cell_cycle_checkpoint

- **Comparison:** (R+45)_vs_(R+1)
- **pathway_rank_metric:** 1.06554657426435
- **average_signed_effect:** 0.4214576180863333
- **overlapping_gene_count:** 3
- **Why it deserves attention:** Exploratory pathway ranking with overlap ≥3 and main-candidate filters passed; interpretation label: moderate positive DNA damage-response signature.
- **Caution:** limited pathway-level signal; small crew sizes possible; not clinical severity.

### OSD-570_PBMC_snRNA | Dendritic Cell | oxidative_stress

- **Comparison:** (R+45)_vs_(R+1)
- **pathway_rank_metric:** 0.6494534676916411
- **average_signed_effect:** -0.3708752685101664
- **overlapping_gene_count:** 3
- **Why it deserves attention:** Exploratory pathway ranking with overlap ≥3 and main-candidate filters passed; interpretation label: weak or unclear signature.
- **Caution:** limited pathway-level signal; small crew sizes possible; not clinical severity.

### OSD-570_PBMC_snRNA | CD16+ Monocyte | DNA_repair

- **Comparison:** (R+45)_vs_(R+1)
- **pathway_rank_metric:** 0.356870675944015
- **average_signed_effect:** 0.154871234092737
- **overlapping_gene_count:** 3
- **Why it deserves attention:** Exploratory pathway ranking with overlap ≥3 and main-candidate filters passed; interpretation label: weak or unclear signature.
- **Caution:** limited pathway-level signal; small crew sizes possible; not clinical severity.

## OSD-569 whole-blood context

OSD-569 is **scoreable** in this pipeline and can show **broad** curated-gene coverage in whole blood, but observed effects are typically **small**, and **no OSD-569 pathway** passes the **main-candidate** threshold in the reference configuration. Treat OSD-569 as **broad whole-blood context**, not the strongest pathway signal on its own.

## OSD-570 PBMC context

OSD-570 tends to carry the **strongest candidate pathway-level signals** in this project when present, but the packaged contrast is **R+45 vs R+1**—**recovery-period timing**, not a direct flight-vs-pre-flight baseline.

## Genes worth closer inspection

Source: `outputs/tables/top_genes_from_pathway_level_signals.csv`. Listed genes are **worth closer inspection** under this scoring framework—they are **not** labeled as damaged genes.

### Gene SOD1 (oxidative_stress)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** CD16+ Monocyte
- **gene_signal:** 2.9531668168657728 | **effect_size:** -0.455083026830571 | **adjusted_p_value:** 3.24121338046725e-07
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene ATR (DNA_repair)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** Other
- **gene_signal:** 2.902442296825451 | **effect_size:** 0.423728977831044 | **adjusted_p_value:** 1.41331456625811e-07
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene RB1 (cell_cycle_checkpoint)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** CD16+ Monocyte
- **gene_signal:** 2.4187760374112472 | **effect_size:** 0.513938035558891 | **adjusted_p_value:** 1.96626784278424e-05
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene TXN (oxidative_stress)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** Dendritic Cell
- **gene_signal:** 1.9483604030749235 | **effect_size:** -0.571782695638896 | **adjusted_p_value:** 0.0003912740028502
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene XRCC5 (DNA_repair)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** Other
- **gene_signal:** 1.5346591973781722 | **effect_size:** 0.302384677445719 | **adjusted_p_value:** 8.41030281282907e-06
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene PRDX1 (oxidative_stress)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** CD16+ Monocyte
- **gene_signal:** 0.7901674765203603 | **effect_size:** -0.320141475413348 | **adjusted_p_value:** 0.003402654435379
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene ATM (DNA_repair)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** CD16+ Monocyte
- **gene_signal:** 0.7778636853818024 | **effect_size:** 0.447717450468995 | **adjusted_p_value:** 0.0183063334724482
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene ATM (cell_cycle_checkpoint)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** CD16+ Monocyte
- **gene_signal:** 0.7778636853818024 | **effect_size:** 0.447717450468995 | **adjusted_p_value:** 0.0183063334724482
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene GPX4 (oxidative_stress)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** CD16+ Monocyte
- **gene_signal:** 0.7529581130578799 | **effect_size:** -0.2965945789497 | **adjusted_p_value:** 0.0028928243574288
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene ATM (DNA_repair)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** Other
- **gene_signal:** 0.4209902990801619 | **effect_size:** 0.514432819027536 | **adjusted_p_value:** 0.151929403420046
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene XRCC6 (DNA_repair)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** CD16+ Monocyte
- **gene_signal:** 0.2927483424502429 | **effect_size:** -0.285821116421898 | **adjusted_p_value:** 0.0945722607581792
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

### Gene XRCC6 (DNA_repair)

- **Dataset:** OSD-570_PBMC_snRNA | **Cell type:** Dendritic Cell
- **gene_signal:** 0.0875409556275018 | **effect_size:** -0.34129255618798 | **adjusted_p_value:** 0.55398973054878
- **Caution:** RNA expression reflects many biology layers (inflammation, stress, immune shifts, composition); does not localize injury or predict cancer.

## Practical use

Frame this project as **triage for astronaut omics follow-up**: it can help decide which **pathways**, **genes**, **immune-cell contexts**, or **assays** deserve closer attention after spaceflight.

Possible follow-up assays include (examples only): **direct sequencing**, **DNA lesion assays**, **comet assay**, **gamma-H2AX staining**, **cytogenetics**, **protein-level validation**, and other orthogonal measures—not inferred from RNA alone.

These signatures may help **prioritize long-term monitoring questions**, but they **do not diagnose cancer** or **predict tissue-specific disease**.

## What not to overclaim

- Does **not** detect DNA damage directly.
- Does **not** identify damaged genes.
- Does **not** predict cancer.
- Does **not** localize injury to a body part.
- RNA expression can reflect **inflammation**, **stress**, **immune shifts**, or **cell-composition** changes—not only DNA repair biology.

---

_Generated by `scripts/build_follow_up_priorities.py` from existing outputs._