# Extended data — supplementary datasets

Summary-level datasets underlying the figures and tables of the manuscript
*Whole-exome sequencing of 12 African populations reveals up to 5.6-fold variation
in predicted drug toxicity risk for essential medicines.*

Raw whole-exome sequencing data are under EGA controlled access
(accession **EGAS00001008456**). The files here are derived, summary-level
outputs from the DeepVariant WES +2 kb-padded pipeline (GRCh38, 127 samples).

## Status of the five promised "Extended data" categories

| # | Dataset | Status | File(s) here |
|---|---------|--------|--------------|
| 1 | Per-population allele frequencies, clinically actionable variants | ⚠️ **needs regeneration** | not yet packaged — see note |
| 2 | The 179 highly differentiated variants (HDVs) | ⚠️ **needs full columns** | not yet packaged — see note |
| 3 | Novel pharmacogene variant catalogue | ⚠️ **definition in flux** | not yet packaged — see note |
| 4 | PyPGx diplotype / metaboliser calls | ✅ included | `pypgx_phenotypes_by_population.csv`, `supp_table_S3_pypgx_phenotypes.csv`, `pgx_star_allele_summary.tsv`, `nat2_phased_acetylator_phenotypes.tsv` |
| 5 | HLA typing results (pharmacorelevant alleles) | ✅ included (summary) | `table9_hla_alleles.csv` |

## Notes on the categories still to package

**1 — Per-population allele frequencies (actionable variants).** The existing
workspace files (`table6_from_genotypes.csv`, `pgx_28_variant_status.csv`,
`table4_pgx_variants_28.csv`) carry **scrambled / off-build coordinates**
(e.g. rs8187710 was listed at chr10:99782821, which is actually rs717620; its
true GRCh38 locus is chr10:99851537). The clean per-population MAFs currently
exist only inside Table 6 of the manuscript `.tex`. Regenerate directly from the
joint VCF before depositing.

**2 — The 179 HDVs.** Only a bare `CHROM:POS:REF:ALT` list exists
(`manuscript_data/hdv/ALL_hdv-…g2_all_all.csv`, 179 rows). A deposit-quality
table needs gene, per-population AF, overall AF, effect class, and the Fisher
exact P added (joinable from `ALL.frq` + the `fisher-test/` outputs on the
cluster). 179 variants span 121 genes (179 records → 173 unique positions → 121
genes).

**3 — Novel pharmacogene variant catalogue.** No clean catalogue file exists yet,
and the count is unsettled: prose says **565**, the regenerated figure caption
says **469**, and an independent strict recompute (dbSNP b156 + gnomAD v4.1
all-pop + AGVP b38, 2026-06-09) gives **470**. The exome-wide novel set was also
recomputed to **16,934** (vs the manuscript's 52,748, which used the looser
gnomAD-AFR-only filter). Settle the definition, then export the catalogue
(key, gene, effect, per-population carriers) from the recompute work dir
`/scratch3/users/mamana/exome_aibst_b38/sweep_checks/novel_v2/`.

## Provenance

- Pipeline: DeepVariant (nf-core/sarek 3.5.1) → GLnexus joint calling → snpEff →
  dbSNP b156 → ClinVar → CADD v1.7 → AlphaMissense → dbNSFP → PyPGx / arcasHLA /
  OptiType. GRCh38, 127 samples.
- Code: this repository (`exome_aibst_b38`).
