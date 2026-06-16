# Extended data — supplementary datasets

Summary-level datasets underlying the figures and tables of the manuscript
*Whole-exome sequencing of 12 African populations reveals up to 4.9-fold variation
in predicted drug toxicity risk for essential medicines.*

Raw whole-exome sequencing data are under EGA controlled access
(accession **EGAS00001008456**). The files here are derived, summary-level
outputs from the DeepVariant WES +2 kb-padded pipeline (GRCh38, 127 samples).

## Status of the five promised "Extended data" categories

| # | Dataset | Status | File(s) here |
|---|---------|--------|--------------|
| 1 | Per-population allele frequencies, clinically actionable variants | ✅ included | `supp_table_actionable_variant_freqs.tsv` |
| 2 | The 173 highly differentiated variants (HDVs) | ✅ included | `supp_table_hdvs.tsv` |
| 3 | Novel pharmacogene variant catalogue | ✅ included | `supp_table_novel_pgx_catalogue.tsv` |
| 4 | PyPGx diplotype / metaboliser calls | ✅ included | `pypgx_phenotypes_by_population.csv`, `supp_table_S3_pypgx_phenotypes.csv`, `pgx_star_allele_summary.tsv`, `nat2_phased_acetylator_phenotypes.tsv` |
| 5 | HLA typing results (pharmacorelevant alleles) | ✅ included (summary) | `table9_hla_alleles.csv` |

## Notes on the packaged tables (1–3, built 2026-06-11)

**1 — `supp_table_actionable_variant_freqs.tsv`** (25 actionable variants, 15
genes). Per-population minor-allele frequencies (%) for the n≥8 populations
(KNL, KNK, NGY, TZA, TZB, SAV, ZWD, ZWS) plus gnomAD AFR/EUR and the AFR-vs-EUR
Fisher exact P, extracted from the manuscript's verified Table 6. NOTE: this
supersedes the buggy `manuscript_data/tables/table6_from_genotypes.csv` /
`pgx_28_variant_status.csv`, which carried **scrambled / off-build coordinates**
(e.g. rs8187710 was wrongly at chr10:99782821 = rs717620; its true GRCh38 locus
is chr10:99851537). `--` = not reliably callable (monomorphic / outside capture
/ insufficient depth), per the manuscript footnotes.

**2 — `supp_table_hdvs.tsv`** (173 highly differentiated variants, 121
genes). rsID, locus, gene, effect, overall AiBST AF, all 12 per-population AFs,
and max pairwise F_ST + the population pair driving it (from
`pgx_max_fst_per_variant.tsv`; 80/173 have an F_ST value, the rest `.`). The 173
variants span 121 unique genes, matching the manuscript's 173 HDVs (deduplicated
from 179 raw records: 6 byte-identical duplicates removed).

**3 — `supp_table_novel_pgx_catalogue.tsv`** (470 novel pharmacogene variants,
all-population definition, ≈ the manuscript's 469). Locus, gene, effect, overall
AF, 12 per-population AFs, carrier populations, and N populations. 86.4% are
single-population (one carrier population). Definition = absent from dbSNP b156 +
gnomAD v4.1 (all populations) + AGVP b38 — the same all-pop definition the
manuscript adopted for Table 2 (17,528 exome-wide / 469 PGx). 410/470 carry
per-population detail (60 lack a per-population freq record, shown as `.`).
Recompute work dir: `/scratch3/users/mamana/exome_aibst_b38/sweep_checks/novel_v2/`.

## Provenance

- Pipeline: DeepVariant (nf-core/sarek 3.5.1) → GLnexus joint calling → snpEff →
  dbSNP b156 → ClinVar → CADD v1.7 → AlphaMissense → dbNSFP → PyPGx / arcasHLA /
  OptiType. GRCh38, 127 samples.
- Code: this repository (`exome_aibst_b38`).
