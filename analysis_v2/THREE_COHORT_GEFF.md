# Current-production effective genome sizes

Status: locally implemented; real-cluster validation pending. No full DA runs.

On 2026-10-01 André verified the historical extraction's mapping and provenance
checksums on the cluster. Mapping:
work/yachida_geff_audit_20260917T144215Z/metaphlan_genome_sizes/metaphlan_genome_sizes.tsv.
SHA256: 649580e792e105fa9144ff8108b93183530a16a26f660b6872afa2fcf899493a.
Its provenance reports mpa_vJan25_CHOCOPhlAnSGB_202503.pkl with database SHA256
e7d23a73a7959b4f41af0bbe403f4b5bbb7c1879528d376d146ee5294515df9a,
58,216 mapped terminal SGBs and PASS. Reuse this mapping with current baselines;
do not reuse the old 507-sample G_eff estimates.

run_three_cohort_geff.sbatch uses one receipt-verified raw baseline per sample,
selected from each cohort's production_seal_v2 and current configured results
root. Expected counts are Yachida 201, Feng 154, Zeller 156. Older results_old_v2
profiles are never searched. Existing upstream outputs and seals are read-only.

Each sample-wide G_eff serves community and independent endpoints. Primary
settings remain terminal SGB rank, exact lineage mapping, 95% eligible-abundance
coverage and zero exclusions. The job rebuilds ten target lengths from frozen
FASTA references and writes per-cohort G_eff tables, input hashes and a final
recursive checksum manifest. Partial failed outputs must not be reused.

Required exports: PROJECT, fresh GEFF_ROOT, ANALYSIS_SIF, GENOME_SIZE_MAPPING
and EXPECTED_MAPPING_SHA256. The mapping must be the previously audited
taxonomy genome-length extraction from the actual production MetaPhlAn
database, not target_genome_sizes.tsv. Establish its database provenance from
the historical extraction records before using its checksum. A checksum alone
proves bytes, not correct database identity. No automatic mapping fallback.

First locate that mapping through the old audit's checksum manifest:

```bash
cat work/three_cohort_geff_audit_20260921T153427Z/geff_primary_cov95/effective_genome_size.sha256
```

Before submission run test_sealed_metaphlan_baselines.py,
test_effective_genome_size.py, test_geff_canonical_integration.py and
test_target_genome_sizes.py using the pinned analysis image. The selectors'
fixture expectations are intentionally one sample, while the production job
always requires full cohort counts. Bash syntax must also pass.

After job completion verify root SHA256SUMS, each cohort SUCCESS/count,
mapping coverage and zero exclusions; compare baseline identities and paths
against the definitive canonical inputs using validate_metaphlan_reference_inputs.py.
This is reference preparation, not a full downstream model run. Inputs outside
/mnt/nfs need an explicit additional Apptainer bind in the site wrapper.
