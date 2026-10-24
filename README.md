# WGS Germline Analysis Pipeline

A containerized germline whole-genome sequencing (WGS) workflow for **single-sample/proband rare-disease analysis**.

The pipeline combines:

- VCF normalization and genotype QC
- HPO-based phenotype representation
- Exomiser phenotype-driven prioritization
- population-frequency filtering
- pathogenicity and functional-impact evidence
- possible compound-heterozygous candidate generation
- secondary-findings candidate generation

> **Important:** This workflow is intended for research, development, benchmarking, and workflow validation. It is not a standalone clinical interpretation system.

---

## Table of Contents

- [Current Scope](#current-scope)
- [Pipeline Overview](#pipeline-overview)
- [Repository Structure](#repository-structure)
- [Docker Image](#docker-image)
- [External Runtime Resources](#external-runtime-resources)
- [Reference Genome](#reference-genome)
- [Input VCF Requirements](#input-vcf-requirements)
- [HPO Input](#hpo-input)
- [External Configuration](#external-configuration)
- [Software and Database Versions](#software-and-database-versions)
- [Build the Docker Image](#build-the-docker-image)
- [Run the Pipeline](#run-the-pipeline)
- [Command-Line Arguments](#command-line-arguments)
- [Pipeline Workflow](#pipeline-workflow)
- [Diagnostic Prioritization](#diagnostic-prioritization)
- [Compound-Heterozygous Candidates](#possible-compound-heterozygous-candidates)
- [Secondary Findings](#secondary-findings)
- [Secondary Input](#secondary-input)
- [Output Structure](#output-structure)
- [Runtime Logging](#runtime-logging)
- [Java Memory](#java-memory)
- [Database Version Manifest](#database-version-manifest)
- [Methods Report](#methods-report)
- [Reproducibility](#reproducibility)
- [Troubleshooting](#troubleshooting)
- [Pre-Submission Checklist](#recommended-pre-submission-checks)
- [Clinical Use Disclaimer](#clinical-use-disclaimer)
- [License](#license)
- [Maintainer](#maintainer)

---

## Current Scope

### Supported

- Germline WGS
- Proband/single-sample analysis
- SNVs
- Small indels
- GRCh38
- Phenotype-driven rare-disease prioritization
- Possible compound-heterozygous candidate generation
- Secondary-findings candidate generation

### Not currently included

- CNV prioritization
- Structural-variant prioritization
- Trio-based inheritance confirmation
- Automatic phasing confirmation
- Final clinical ACMG/AMP interpretation

---

## Pipeline Overview

```text
Input VCF
   |
   v
Primary chromosome filtering
   |
   v
Variant normalization + multiallelic splitting
   |
   v
PASS filtering
   |
   v
Remove homozygous-reference variants
   |
   v
DP / GQ / allele-balance QC
   |
   +----------------------------+
   |                            |
   v                            v
Generate Exomiser YAML      QC VCF extraction
from HPO terms              to TSV
   |
   v
Run Exomiser 15.1.0
   |
   v
Extract Exomiser TSV + JSONL evidence
   |
   v
Merge Exomiser TSV + JSON evidence
   |
   v
Merge VCF + Exomiser evidence
   |
   v
Final prioritization
   |
   +-----------------------------+
   |              |              |
   v              v              v
Top-N         Possible       Secondary
variants      compound-het   findings
              pairs          candidates
```

---

## Repository Structure

Recommended layout:

```text id="qj0wq8"
wgs_pipeline/
├── Dockerfile
├── run_wgs_pipeline.sh
├── README.md
│
├── scripts/
│   ├── yml_hpo_parser.py
│   ├── wgs_vcf_extract.py
│   ├── exomiser_json_extract.py
│   ├── exomiser_json_tsv_merge.py
│   ├── wgs_vcf_exomiser_merge.py
│   └── wgs_variants_prioritizer.py
│
└── config-example/
    ├── wgs_pipeline.env.example
    ├── hg38application.properties.example
    ├── exomiser_wgs_template.yml
    ├── secondary_findings_master.tsv
    └── database_versions.tsv
```

---

## Docker Image

The Docker image contains:

- Ubuntu 24.04
- OpenJDK 21
- Python 3
- Polars
- bcftools
- tabix
- Exomiser CLI 15.1.0
- Pipeline Python scripts
- Pipeline runner script

User-editable references, databases, and configuration files are intentionally kept **outside the Docker image**.

---

## External Runtime Resources

The following are mounted at runtime:

- Input VCF
- Output directory
- Reference FASTA
- Exomiser data
- Exomiser properties file
- Exomiser YAML template
- Secondary-findings master file
- Pipeline environment configuration

This design allows reference/database changes without rebuilding the Docker image.

---

## Reference Genome

The workflow uses:

```text
GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

> **Reference compatibility:** Chromosome/contig names in the reference FASTA must match those used in the input VCF (for example, `1` vs `chr1`, and `MT` vs `chrM`). Ensure the naming convention is consistent before running the pipeline.

---

## Input VCF Requirements

The input VCF must:

- be bgzip-compressed
- end in `.vcf.gz`
- contain one proband/sample
- contain the FORMAT fields:
  - `GT`
  - `DP`
  - `GQ`
  - `AD`

Example:

```text
sample.vcf.gz
```

---

## HPO Input

HPO identifiers are supplied as a comma-separated string.

Example:

```text
HP:0001250,HP:0001263
```

These are inserted into a case-specific Exomiser YAML file using:

```text
scripts/yml_hpo_parser.py
```

---

## External Configuration

The workflow uses one external configuration file:

```text
wgs_pipeline.env
```

Example:

```bash
REFERENCE_FASTA="/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

EXOMISER_PROPERTIES="/config/hg38application.properties"

EXOMISER_TEMPLATE="/config/exomiser_wgs_template.yml"

SECONDARY_GENE_LIST="/config/secondary_findings_master.tsv"

EXOMISER_DATA_DIR="/exomiser-data"

THREADS=8

MIN_DP=10
MIN_GQ=20

HET_AB_MIN=0.20
HET_AB_MAX=0.80

TOP_N=25

PRIMARY_CHROMS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"

EXOMISER_ASSEMBLY="hg38"

JAVA_XMS="4g"
JAVA_XMX="16g"
```

---

## Software and Database Versions

| Resource | Version / release | Purpose |
|---|---|---|
| Reference genome | `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` | GRCh38 normalization/reference |
| Exomiser | `15.1.0` | Phenotype-driven prioritization |
| Exomiser data | `2602` | Phenotype, disease, population, inheritance, and variant evidence |
| CADD | `v1.7` | Pathogenicity prediction |
| ClinGen | `3.3` | Secondary-finding and gene-disease guidance |
| HPO | Included in Exomiser `2602` | Phenotype representation |
| gnomAD | Included in Exomiser `2602` | Population frequency evidence |
| ClinVar | Included in Exomiser `2602` | Clinical significance evidence |
| REVEL | Included through Exomiser `2602` where available | Missense prediction |
| SpliceAI | Included through Exomiser `2602` where available | Splice-effect prediction |
| Java | OpenJDK 21 | Exomiser runtime |

The machine-readable manifest is maintained in:

```text
config-example/database_versions.tsv
```

---

## Build the Docker Image

From the repository root:

```bash
docker build -t vgen/wgs-pipeline:0.1 .
```

---

## Run the Pipeline

Example:

```bash
docker run --rm -it \
    -v "$PWD/input:/input:ro" \
    -v "$PWD/output:/output" \
    -v "$PWD/config:/config:ro" \
    -v "$PWD/reference:/reference:ro" \
    -v "$PWD/exomiser-data:/exomiser-data:ro" \
    vgen/wgs-pipeline:0.1 \
    --config /config/wgs_pipeline.env \
    --input-vcf /input/sample.vcf.gz \
    --case-id CASE001 \
    --hpo-ids "HP:0001250,HP:0001263" \
    --output-dir /output
```

---

## Command-Line Arguments

| Argument | Required | Description |
|---|---:|---|
| `--config` | Yes | External pipeline configuration file |
| `--input-vcf` | Yes | Input bgzip-compressed VCF |
| `--case-id` | Yes | Case/sample identifier |
| `--hpo-ids` | Yes | Comma-separated HPO IDs |
| `--output-dir` | Yes | Output directory |
| `--top N` | No | Overrides `TOP_N` from configuration |

---

## Pipeline Workflow

### Step 1: Primary Chromosome Filtering

Default chromosomes:

```text
1-22, X, Y, MT
```

Configured through:

```text
PRIMARY_CHROMS
```

### Step 2: Variant Normalization

Variants are normalized against the configured GRCh38 FASTA.

Conceptually:

```bash
bcftools norm \
    -f "$REFERENCE_FASTA" \
    -m -any
```

This performs:

- reference-based normalization
- multiallelic splitting

### Step 3: PASS Filtering

Only:

```text
FILTER=PASS
```

variants are retained.

PASS filtering is performed **after normalization**.

### Step 4: Homozygous-Reference Removal

Homozygous-reference genotypes are removed before tertiary analysis.

### Step 5: Genotype QC

Default thresholds:

```text
DP >= 10
GQ >= 20
```

For heterozygous variants:

```text
AB = ALT_AD / (REF_AD + ALT_AD)
```

Accepted range:

```text
0.20 <= AB <= 0.80
```

Homozygous-alternate and haploid-alternate variants are evaluated using DP/GQ.

### Step 6: Generate Case-Specific Exomiser YAML

The pipeline uses:

```text
scripts/yml_hpo_parser.py
```

to insert supplied HPO IDs into the configured Exomiser YAML template.

### Step 7: Run Exomiser

The workflow uses:

```text
Exomiser CLI 15.1.0
Exomiser data release 2602
```

Exomiser runs directly inside the Docker container.

No nested Docker is used.

### Step 8: QC VCF Extraction

`wgs_vcf_extract.py` extracts:

```text
CHROM
POS
REF
ALT
QUAL
FILTER
INFO
GT
ZYGOSITY
AD
VAF
DP
GQ
PL
```

### Step 9: Exomiser JSONL Extraction

`exomiser_json_extract.py` extracts available evidence including:

- population frequencies
- pathogenicity scores
- phenotype scores
- inheritance evidence
- contributing-variant status
- whitelist status
- ClinVar evidence
- Exomiser ACMG-related fields

### Step 10: Exomiser Evidence Merge

`exomiser_json_tsv_merge.py` merges:

```text
Exomiser variants.tsv
+
Exomiser JSONL-derived TSV
```

using:

```text
CONTIG + START + REF + ALT
```

### Step 11: VCF + Exomiser Merge

`wgs_vcf_exomiser_merge.py` joins:

```text
VCF:
CHROM + POS + REF + ALT
```

against:

```text
Exomiser:
CONTIG + START + REF + ALT
```

Mapping:

```text
VCF.CHROM = Exomiser.CONTIG
VCF.POS   = Exomiser.START
VCF.REF   = Exomiser.REF
VCF.ALT   = Exomiser.ALT
```

Chromosome names are normalized before joining.

Examples:

```text
chr1 -> 1
chrX -> X
chrY -> Y
chrM -> MT
M    -> MT
```

The final merged table contains canonical downstream fields:

```text
CONTIG
START
END
REF
ALT
```

and retains:

```text
CHROM
POS
```

for VCF provenance.

The merge also writes:

```text
EXOMISER_MATCH_STATUS
```

with values:

```text
MATCHED
NOT_MATCHED
```

### Step 12: Final Variant Prioritization

`wgs_variants_prioritizer.py` generates:

- Top-N diagnostic candidates
- Possible compound-heterozygous pairs
- Secondary-findings candidates

---

## Diagnostic Prioritization

### Population Frequency

The workflow calculates the maximum configured gnomAD subpopulation allele frequency.

Default thresholds:

| Inheritance | Maximum AF |
|---|---:|
| Autosomal recessive | `< 0.01` |
| Autosomal dominant | `< 0.001` |
| Unknown / other | `< 0.01` |

Variants with all configured gnomAD frequency fields missing are allowed to continue.

### High-Impact Consequences

Current high-impact consequences:

```text
stop_gained
frameshift_variant
start_lost
stop_lost
splice_acceptor_variant
splice_donor_variant
```

### Predictor Evidence

The prioritizer can use:

- CADD v1.7
- REVEL
- SpliceAI

Predictor evidence is used for candidate prioritization and is not treated as standalone clinical classification evidence.

### Rescue Logic

Variants may be retained through strong existing evidence such as:

- whitelist status
- ClinVar pathogenic/likely pathogenic evidence
- Exomiser pathogenic/likely pathogenic evidence

### Final Candidate Logic

Conceptually:

```text
POPULATION_PASS
AND
EXOMISER_RELEVANT
AND
(
    RESCUE
    OR HIGH_IMPACT
    OR PREDICTOR_PASS
)
```

### Top-N Ranking

Variants are ranked using:

1. Exomiser rank
2. Whitelist status
3. Contributing-variant status
4. Exomiser gene combined score
5. Exomiser variant score
6. Exomiser phenotype score
7. Population frequency

Variants are deduplicated by:

```text
CONTIG + START + REF + ALT
```

---

## Possible Compound-Heterozygous Candidates

For autosomal-recessive disease models, qualifying heterozygous variants are grouped by:

- gene
- disease context
- inheritance mode

When at least two qualifying distinct variants exist, candidate pairs are generated.

Output labels:

```text
PAIR_STATUS = POSSIBLE_COMPOUND_HET
PHASE_STATUS = UNKNOWN
```

Because this is proband-only analysis, the workflow does not claim that two variants are confirmed to be in trans.

---

## Secondary Findings

Secondary findings are evaluated separately from phenotype-driven diagnostic prioritization.

The secondary-findings branch does **not** require:

- high Exomiser rank
- HPO similarity
- `CONTRIBUTING_VARIANT` status
- inclusion in the Top-N diagnostic list

The workflow uses:

```text
ClinGen SF v3.3
```

together with the configured secondary-findings master file.

Candidate states can include:

```text
SF_CANDIDATE_REVIEW
PHASING_REQUIRED
SF_MANUAL_REVIEW
```

Secondary findings remain candidate-level outputs and require appropriate expert review before clinical reporting.

---

## Secondary Input

The prioritizer supports the optional argument:

```text
--secondary-input
```

The preferred secondary input is a **complete QC-passed annotated variant table** independent of Exomiser phenotype filtering.

Minimum useful fields:

```text
GENE_SYMBOL
CONTIG
START
REF
ALT
GENOTYPE
CLINVAR_PRIMARY_INTERPRETATION
```

If no separate full annotated table is available, the pipeline can fall back to the merged Exomiser/VCF table.

This fallback may be less complete for secondary-findings analysis because Exomiser is phenotype-driven.

---

## Output Structure

Typical output:

```text
output/
├── work/
│   ├── CASE_primary_chr.vcf.gz
│   ├── CASE_normalized_split.vcf.gz
│   ├── CASE_pass.vcf.gz
│   ├── CASE_nonref.vcf.gz
│   └── CASE_exomiser.yml
│
├── results/
│   ├── CASE_qc_filtered.vcf.gz
│   ├── CASE_vcf.tsv
│   ├── CASE_exomiser.variants.tsv
│   ├── CASE_exomiser.jsonl
│   ├── CASE_exomiser_details.tsv
│   ├── CASE_exomiser_enriched.tsv
│   └── CASE_wgs_merged.tsv
│
├── prioritized/
│   ├── CASE_top25.tsv
│   ├── CASE_possible_compound_het.tsv
│   └── CASE_secondary_findings.tsv
│
└── logs/
    └── CASE_pipeline.log
```

---

## Runtime Logging

The pipeline records timestamps for analysis stages.

The runner can also record:

- pipeline start time
- pipeline end time
- total elapsed seconds
- total wall-clock runtime

Example:

```text
Runtime:
  Start time     : 2026-10-24 09:01:14
  End time       : 2026-10-24 09:47:36
  Total seconds  : 2782
  Total runtime  : 0h 46m 22s
```

Runtime depends on:

- WGS variant count
- Exomiser data
- CPU allocation
- available memory
- storage performance

---

## Java Memory

The Java heap is configured using separate values:

```bash
JAVA_XMS="4g"
JAVA_XMX="16g"
```

Java should be invoked as:

```bash
java \
    "-Xms${JAVA_XMS}" \
    "-Xmx${JAVA_XMX}" \
    ...
```

This prevents both JVM options from being interpreted as one argument.

---

## Database Version Manifest

The repository should include:

```text
config-example/database_versions.tsv
```

It records:

- resource
- version/release
- source
- purpose
- notes

Current major resources:

```text
Reference     GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
Exomiser      15.1.0
Exomiser data 2602
CADD          v1.7
ClinGen       3.3
```

---
## Methods Report

The Track 1 methods report is submitted separately as part of the hackathon submission.

It includes:

- detailed method description
- public data resources
- AI-assistance disclosure
- compound-heterozygous handling
- secondary-findings handling
- runtime and cost information
- method abstract
- strengths and limitations

---

## Reproducibility

To reproduce the workflow:

1. Clone the repository.
2. Obtain the required external reference and Exomiser data.
3. Create/edit `wgs_pipeline.env`.
4. Build the Docker image.
5. Mount input, output, config, reference, and Exomiser-data directories.
6. Run the pipeline.
7. Review logs and final outputs.

---

## Troubleshooting

### 1. Reference mismatch

If `bcftools norm` reports reference mismatches, verify:

- VCF assembly
- FASTA assembly
- chromosome naming
- reference accession/build

### 2. No Exomiser matches after merge

Compare:

```text
VCF:
CHROM
POS
REF
ALT
```

against:

```text
Exomiser:
CONTIG
START
REF
ALT
```

Check:

- chromosome naming
- coordinate system
- REF/ALT representation
- normalization status

### 3. Missing `CONTIG` / `START` downstream

The final merged table must contain:

```text
CONTIG
START
END
REF
ALT
```

The current merge script explicitly creates these canonical columns.

### 4. Missing `AD`

The heterozygous allele-balance filter requires:

```text
FORMAT/AD
```

### 5. Exomiser does not start

Verify:

- Java version
- Exomiser JAR
- properties file
- Exomiser data path
- JVM heap parameters

### 6. Exomiser `FileNotFoundException` for linked files

If symbolic links are used inside mounted directories, the symlink target must also be visible inside the Docker container.

Check:

```bash
readlink -f /exomiser-data/<file>
```

### 7. No secondary findings

Review:

- secondary-findings gene overlap
- `GENE_SYMBOL` consistency
- ClinVar P/LP fields
- Exomiser P/LP fields
- inheritance state
- whether a complete `--secondary-input` table was supplied

---

## Recommended Pre-Submission Checks

- [ ] Docker image builds from a clean checkout
- [ ] All Python scripts are committed
- [ ] `run_wgs_pipeline.sh` is committed
- [ ] Example configuration is committed
- [ ] `database_versions.tsv` is committed
- [ ] `track1_report.md` is committed
- [ ] No patient identifiers are present
- [ ] No API keys or credentials are present
- [ ] Restricted/licensed databases are not committed improperly
- [ ] Runtime is recorded for at least one complete case
- [ ] Final software versions are recorded
- [ ] Final Docker image tag is documented

---

## Clinical Use Disclaimer

This workflow is intended for:

- research
- development
- benchmarking
- validation
- bioinformatics workflow implementation

Any clinical use requires appropriate:

- analytical validation
- quality management
- expert variant review
- disease-specific interpretation
- segregation/phasing where relevant
- orthogonal confirmation where required
- clinical reporting procedures

Automated prioritization output must not be interpreted as a standalone clinical diagnosis.

---

## License


This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Maintainer

**Vgenomics Pvt. Ltd.**
