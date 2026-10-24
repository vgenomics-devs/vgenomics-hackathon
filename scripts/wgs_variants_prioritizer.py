#!/usr/bin/env python3
import argparse
import itertools
import sys
from pathlib import Path
import polars as pl

# ============================================================
# Configuration
# ============================================================
GNOMAD_AF_COLUMNS = [
    "GNOMAD_E_AFR_AF",
    "GNOMAD_E_AMR_AF",
    "GNOMAD_E_EAS_AF",
    "GNOMAD_E_NFE_AF",
    "GNOMAD_E_SAS_AF",
    "GNOMAD_G_AFR_AF",
    "GNOMAD_G_AMR_AF",
    "GNOMAD_G_EAS_AF",
    "GNOMAD_G_NFE_AF",
    "GNOMAD_G_SAS_AF",
]

VARIANT_KEY = [
    "CONTIG",
    "START",
    "REF",
    "ALT",
]

# ============================================================
# Final Top-N output cleanup
#
# These columns are still used internally where required, but
# are removed from the final Top-N TSV because they are either
# duplicate coordinates, redundant/raw fields, or internal
# prioritization/debug flags.
# ============================================================
TOP_OUTPUT_DROP_COLUMNS = [
    "#RANK",
    "CONTIG",
    "START",
    "END",
    "EXOMISER_MATCH_STATUS",
    "ID",
    "VCF_ID",
    "CHANGE_LENGTH",
    "EXOMISER_QUAL",
    "EXOMISER_FILTER",
    "GENOTYPE",
    "EXOMISER_JSON_GENES",
    "POPULATION_PASS",
    "RESCUE_VARIANT",
    "HIGH_IMPACT",
    "PREDICTOR_PASS",
    "EXOMISER_RELEVANT",
    "KEEP_VARIANT",
    "INFO",
    "PL",
    "P-VALUE",
    "MAX_FREQ_SOURCE",
    "MAX_FREQ",
    "ALL_FREQ",
    "MAX_PATH_SOURCE",
    "MAX_PATH",
    "ALL_PATH",
    "PASSED_FILTER_TYPES",
    "COMPATIBLE_INHERITANCE_MODES",
    "CONTRIBUTING_INHERITANCE_MODES",
    "FINAL_ELIGIBLE",
    "_AR_HET_VARIANT_COUNT",
    "_SCORE_CLINVAR_PLP",
    "_SCORE_EXOMISER_PLP",
    "_SCORE_CLINVAR_CONFLICT",
    "_SCORE_WEAK_NONCODING",
    "_SCORE_SPLICE_SUPPORT",
    "_SCORE_FREQUENCY_BONUS",
    "_SCORE_INHERITANCE_COMPONENT",
    "_EXOMISER_RANK_BONUS",
]

# Preferred leading columns in the final Top-N output.
# Any other retained evidence columns follow in their original order.
TOP_OUTPUT_FIRST_COLUMNS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "INTERNAL_PRIORITY_SCORE",
    "GENE_SYMBOL",
    "ENTREZ_GENE_ID",
    "MOI",
    "RS_ID",
    "GT",
    "ZYGOSITY",
    "AD",
    "VAF",
    "DP",
    "GQ",
    "QUAL",
    "FILTER",
    "FUNCTIONAL_CLASS",
    "HGVS",
    "EXOMISER_GENE_COMBINED_SCORE",
    "EXOMISER_GENE_PHENO_SCORE",
    "EXOMISER_GENE_VARIANT_SCORE",
    "EXOMISER_VARIANT_SCORE",
    "CONTRIBUTING_VARIANT",
    "WHITELIST_VARIANT",
    "EXOMISER_ACMG_CLASSIFICATION",
    "EXOMISER_ACMG_EVIDENCE",
    "EXOMISER_ACMG_DISEASE_ID",
    "EXOMISER_ACMG_DISEASE_NAME",
    "CLINVAR_VARIATION_ID",
    "CLINVAR_PRIMARY_INTERPRETATION",
    "CLINVAR_STAR_RATING",
    "MAX_GNOMAD_SUBPOP_AF",
    "CADD_RAW_SCORE",
    "CADD_SCALED_SCORE",
    "REVEL_SCORE",
    "SPLICE_AI_SCORE",
    "FREQUENCY_SCORE",
    "PATHOGENICITY_SCORE",
    "PREDICTED_PATHOGENIC",
]

HIGH_IMPACT_CLASSES = [
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "stop_lost",
    "splice_acceptor_variant",
    "splice_donor_variant",
]


WEAK_NONCODING_CLASSES = [
    "intergenic_variant",
    "coding_transcript_intron_variant",
    "intron_variant",
    "5_prime_utr_intron_variant",
    "3_prime_utr_intron_variant",
    "5_prime_utr_exon_variant",
    "3_prime_utr_exon_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "regulatory_region_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
]


CLINVAR_RESCUE = [
    "PATHOGENIC",
    "LIKELY_PATHOGENIC",
    "PATHOGENIC_LIKELY_PATHOGENIC",
    "PATHOGENIC_OR_LIKELY_PATHOGENIC",
]

EXOMISER_RESCUE = [
    "PATHOGENIC",
    "LIKELY_PATHOGENIC",
    "PATHOGENIC_LIKELY_PATHOGENIC",
]

HET_GT = [
    "0/1",
    "1/0",
    "0|1",
    "1|0",
]

HOM_ALT_GT = [
    "1/1",
    "1|1",
]

HEMIZYGOUS_ALT_GT = [
    "1",
]

AD_MOI = [
    "AD",
    "AUTOSOMAL DOMINANT",
]

AR_MOI = [
    "AR",
    "AUTOSOMAL RECESSIVE",
]

XL_MOI = [
    "XL",
    "X-LINKED",
    "X LINKED",
    "XLR",
    "XLD",
]

# ============================================================
# Arguments
# ============================================================

def parse_arguments():

    parser = argparse.ArgumentParser(
        description=(
            "WGS germline variant prioritizer. Produces Top-N "
            "diagnostic variants, possible compound heterozygous "
            "pairs, and ACMG secondary-finding candidates."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Merged Exomiser + VCF TSV used for diagnostic prioritization."
    )

    parser.add_argument(
        "--secondary-input",
        required=False,
        default=None,
        help=(
            "Optional full QC-passed annotated WGS TSV for secondary "
            "findings. If omitted, --input is used."
        )
    )

    parser.add_argument(
        "--secondary-gene-list",
        required=True,
        help=(
            "Merged ACMG/ClinGen secondary findings master TSV. "
            "Expected columns: GENE_SYMBOL, INHERITANCE, "
            "REPORTABLE_AS_SF, VARIANTS_TO_REPORT."
        )
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--prefix",
        required=True,
        help="Case ID / output prefix."
    )

    parser.add_argument(
        "--top",
        type=int,
        default=25,
        help="Number of prioritized unique variants. Default: 25."
    )

    return parser.parse_args()

# ============================================================
# Utilities
# ============================================================

def validate_columns(df, required, label):

    missing = [
        col
        for col in required
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{label} missing required column(s): "
            + ", ".join(missing)
        )

def normalize_string(column):

    return (
        pl.col(column)
        .cast(pl.String)
        .str.strip_chars()
    )

def normalize_classification(column):

    return (
        pl.col(column)
        .cast(pl.String)
        .fill_null("")
        .str.strip_chars()
        .str.to_uppercase()
        .str.replace_all(r"[^A-Z0-9]+", "_")
        .str.strip_chars("_")
    )

def plp_expr(column):

    value = (
        pl.col(column)
        .fill_null("")
        .cast(pl.String)
        .str.to_uppercase()
    )

    return (
        value.str.contains(
            r"PATHOGENIC|LIKELY[_ /-]*PATHOGENIC"
        )
        &
        ~value.str.contains(
            r"BENIGN|UNCERTAIN|CONFLICT"
        )
    )

def read_tsv(path):

    # Read all TSV columns initially as strings.
    #
    # This prevents Polars from inferring chromosome/contig columns
    # such as CONTIG or CHROM as Int64 when early rows contain only
    # autosomes (1,2,3,...) and later failing on X, Y, or MT.
    #
    # Numeric fields are explicitly cast later in prepare_dataframe()
    # and prepare_secondary_input().
    return pl.read_csv(
        path,
        separator="\t",
        infer_schema=False,
        null_values=[
            ".",
            "",
            "NA",
            "N/A",
        ],
    )

# ============================================================
# Main diagnostic dataframe preparation
# ============================================================

def prepare_dataframe(df):

    for column in GNOMAD_AF_COLUMNS:

        if column not in df.columns:

            df = df.with_columns(
                pl.lit(None)
                .cast(pl.Float64)
                .alias(column)
            )

    numeric_columns = [
        "#RANK",
        "START",
        "CONTRIBUTING_VARIANT",
        "WHITELIST_VARIANT",
        "EXOMISER_GENE_COMBINED_SCORE",
        "EXOMISER_GENE_PHENO_SCORE",
        "EXOMISER_GENE_VARIANT_SCORE",
        "EXOMISER_VARIANT_SCORE",
        "CADD_RAW_SCORE",
        "REVEL_SCORE",
        "SPLICE_AI_SCORE",
        *GNOMAD_AF_COLUMNS,
    ]

    casts = []

    for column in numeric_columns:

        if column in df.columns:

            casts.append(
                pl.col(column)
                .cast(pl.Float64, strict=False)
                .alias(column)
            )

    if casts:
        df = df.with_columns(casts)

    string_columns = [
        "MOI",
        "GENE_SYMBOL",
        "FUNCTIONAL_CLASS",
        "CONTIG",
        "REF",
        "ALT",
        "GENOTYPE",
    ]

    strings = []

    for column in string_columns:

        if column in df.columns:

            strings.append(
                normalize_string(column)
                .alias(column)
            )

    if strings:
        df = df.with_columns(strings)

    df = df.with_columns(

        pl.col("MOI")
        .str.to_uppercase(),

        pl.col("GENE_SYMBOL")
        .str.to_uppercase(),

        pl.col("FUNCTIONAL_CLASS")
        .str.to_lowercase(),

        pl.col("REF")
        .str.to_uppercase(),

        pl.col("ALT")
        .str.to_uppercase(),

        normalize_classification(
            "CLINVAR_PRIMARY_INTERPRETATION"
        ).alias(
            "CLINVAR_PRIMARY_INTERPRETATION"
        ),

        normalize_classification(
            "EXOMISER_ACMG_CLASSIFICATION"
        ).alias(
            "EXOMISER_ACMG_CLASSIFICATION"
        ),
    )

    return df

# ============================================================
# Secondary input preparation
# ============================================================

def prepare_secondary_input(df):

    required = [
        "GENE_SYMBOL",
        "CONTIG",
        "START",
        "REF",
        "ALT",
        "GENOTYPE",
        "CLINVAR_PRIMARY_INTERPRETATION",
    ]

    validate_columns(
        df,
        required,
        "Secondary findings variant input"
    )

    string_columns = [
        "GENE_SYMBOL",
        "CONTIG",
        "REF",
        "ALT",
        "GENOTYPE",
    ]

    expressions = []

    for column in string_columns:

        expressions.append(
            pl.col(column)
            .cast(pl.String)
            .str.strip_chars()
            .alias(column)
        )

    df = df.with_columns(expressions)

    df = df.with_columns(

        pl.col("GENE_SYMBOL")
        .str.to_uppercase(),

        pl.col("REF")
        .str.to_uppercase(),

        pl.col("ALT")
        .str.to_uppercase(),

        pl.col("START")
        .cast(pl.Int64, strict=False),

        normalize_classification(
            "CLINVAR_PRIMARY_INTERPRETATION"
        ).alias(
            "CLINVAR_PRIMARY_INTERPRETATION"
        ),
    )

    if "EXOMISER_ACMG_CLASSIFICATION" not in df.columns:

        df = df.with_columns(
            pl.lit(None)
            .cast(pl.String)
            .alias("EXOMISER_ACMG_CLASSIFICATION")
        )

    else:

        df = df.with_columns(

            normalize_classification(
                "EXOMISER_ACMG_CLASSIFICATION"
            ).alias(
                "EXOMISER_ACMG_CLASSIFICATION"
            )
        )

    return df

# ============================================================
# Diagnostic prioritization flags
# ============================================================

def add_prioritization_flags(df):

    df = df.with_columns(

        pl.max_horizontal(
            [
                pl.col(column)
                for column in GNOMAD_AF_COLUMNS
            ]
        )
        .alias("MAX_GNOMAD_SUBPOP_AF")
    )

    df = df.with_columns(

        pl.when(
            pl.col("MAX_GNOMAD_SUBPOP_AF").is_null()
        )
        .then(True)

        .when(
            pl.col("MOI") == "AR"
        )
        .then(
            pl.col("MAX_GNOMAD_SUBPOP_AF") < 0.01
        )

        .when(
            pl.col("MOI") == "AD"
        )
        .then(
            pl.col("MAX_GNOMAD_SUBPOP_AF") < 0.001
        )

        .otherwise(
            pl.col("MAX_GNOMAD_SUBPOP_AF") < 0.01
        )
        .alias("POPULATION_PASS")
    )

    df = df.with_columns(

        (
            (pl.col("WHITELIST_VARIANT") == 1)
            |
            pl.col(
                "CLINVAR_PRIMARY_INTERPRETATION"
            ).is_in(CLINVAR_RESCUE)
            |
            pl.col(
                "EXOMISER_ACMG_CLASSIFICATION"
            ).is_in(EXOMISER_RESCUE)
        )
        .alias("RESCUE_VARIANT")
    )

    df = df.with_columns(

        pl.col("FUNCTIONAL_CLASS")
        .is_in(HIGH_IMPACT_CLASSES)
        .alias("HIGH_IMPACT")
    )

    # NOTE:
    # CADD_RAW_SCORE >= 15 is retained from your existing logic.
    # Verify whether this column actually represents PHRED-like CADD.
    df = df.with_columns(

        (
            (
                pl.col("CADD_RAW_SCORE")
                .fill_null(float("-inf"))
                >= 15
            )
            |
            (
                pl.col("REVEL_SCORE")
                .fill_null(float("-inf"))
                >= 0.644
            )
            |
            (
                pl.col("SPLICE_AI_SCORE")
                .fill_null(float("-inf"))
                > 0.1
            )
        )
        .alias("PREDICTOR_PASS")
    )

    df = df.with_columns(

        (
            (pl.col("CONTRIBUTING_VARIANT") == 1)
            |
            (pl.col("WHITELIST_VARIANT") == 1)
            |
            pl.col("RESCUE_VARIANT")
        )
        .alias("EXOMISER_RELEVANT")
    )

    df = df.with_columns(

        (
            pl.col("POPULATION_PASS")
            &
            pl.col("EXOMISER_RELEVANT")
            &
            (
                pl.col("RESCUE_VARIANT")
                |
                pl.col("HIGH_IMPACT")
                |
                pl.col("PREDICTOR_PASS")
            )
        )
        .alias("KEEP_VARIANT")
    )

    return df


# ============================================================
# Internal priority score
# ============================================================

def add_internal_priority_score(df):
    """
    Add a general, case-agnostic INTERNAL_PRIORITY_SCORE and
    FINAL_ELIGIBLE flag.

    IMPORTANT:
      - INTERNAL_PRIORITY_SCORE is NOT a probability.
      - It is NOT an ACMG/AMP classification.
      - No gene-specific, disease-specific, or case-specific rules are used.
      - Exomiser #RANK is used internally but is removed from the final TSV.

    The layer has two stages:

      1. FINAL_ELIGIBLE gate
         Removes weak background candidates before final ranking.

      2. INTERNAL_PRIORITY_SCORE
         Ranks the remaining candidates using phenotype, Exomiser evidence,
         inheritance, consequence, population rarity, pathogenic evidence,
         predictors, and a bounded Exomiser-rank bonus.
    """

    # --------------------------------------------------------
    # Count distinct heterozygous variants per gene for AR models.
    #
    # This provides possible biallelic support only.
    # It does not establish trans phase.
    # --------------------------------------------------------
    ar_het_counts = (
        df
        .filter(
            pl.col("MOI").is_in(AR_MOI)
            &
            pl.col("GENOTYPE").is_in(HET_GT)
            &
            pl.col("POPULATION_PASS")
        )
        .group_by("GENE_SYMBOL")
        .agg(
            pl.struct(VARIANT_KEY)
            .n_unique()
            .alias("_AR_HET_VARIANT_COUNT")
        )
    )

    df = df.join(
        ar_het_counts,
        on="GENE_SYMBOL",
        how="left",
    )

    df = df.with_columns(
        pl.col("_AR_HET_VARIANT_COUNT")
        .fill_null(0)
        .cast(pl.Int64)
    )

    # --------------------------------------------------------
    # Evidence helper flags.
    # --------------------------------------------------------
    df = df.with_columns(
        plp_expr(
            "CLINVAR_PRIMARY_INTERPRETATION"
        ).alias("_SCORE_CLINVAR_PLP"),

        plp_expr(
            "EXOMISER_ACMG_CLASSIFICATION"
        ).alias("_SCORE_EXOMISER_PLP"),

        (
            pl.col("CLINVAR_PRIMARY_INTERPRETATION")
            .fill_null("")
            .cast(pl.String)
            .str.to_uppercase()
            .str.contains(r"BENIGN|CONFLICT")
        ).alias("_SCORE_CLINVAR_CONFLICT"),

        pl.col("FUNCTIONAL_CLASS")
        .fill_null("")
        .is_in(WEAK_NONCODING_CLASSES)
        .alias("_SCORE_WEAK_NONCODING"),

        (
            pl.col("SPLICE_AI_SCORE")
            .fill_null(float("-inf"))
            > 0.1
        ).alias("_SCORE_SPLICE_SUPPORT"),
    )

    # --------------------------------------------------------
    # General rarity component.
    #
    # Missing AF is treated as unknown/rare-supporting but does
    # not override the final eligibility gate.
    # --------------------------------------------------------
    df = df.with_columns(
        pl.when(
            pl.col("MAX_GNOMAD_SUBPOP_AF").is_null()
        )
        .then(4.0)

        .when(
            pl.col("MAX_GNOMAD_SUBPOP_AF") < 0.0001
        )
        .then(4.0)

        .when(
            pl.col("MAX_GNOMAD_SUBPOP_AF") < 0.001
        )
        .then(3.0)

        .when(
            pl.col("MAX_GNOMAD_SUBPOP_AF") < 0.01
        )
        .then(1.0)

        .otherwise(0.0)

        .alias("_SCORE_FREQUENCY_BONUS")
    )

    # --------------------------------------------------------
    # General inheritance component.
    #
    # Genotype alone cannot rescue a weak variant.
    # This component is applied only after eligibility is assessed.
    # --------------------------------------------------------
    df = df.with_columns(
        pl.when(
            pl.col("MOI").is_in(AD_MOI)
            &
            (
                pl.col("GENOTYPE").is_in(HET_GT)
                |
                pl.col("GENOTYPE").is_in(HOM_ALT_GT)
            )
        )
        .then(4.0)

        .when(
            pl.col("MOI").is_in(AR_MOI)
            &
            pl.col("GENOTYPE").is_in(HOM_ALT_GT)
        )
        .then(7.0)

        .when(
            pl.col("MOI").is_in(AR_MOI)
            &
            pl.col("GENOTYPE").is_in(HET_GT)
            &
            (pl.col("_AR_HET_VARIANT_COUNT") >= 2)
        )
        .then(6.0)

        .when(
            pl.col("MOI").is_in(AR_MOI)
            &
            pl.col("GENOTYPE").is_in(HET_GT)
            &
            (pl.col("_AR_HET_VARIANT_COUNT") < 2)
        )
        .then(-6.0)

        .when(
            pl.col("MOI").is_in(XL_MOI)
            &
            (
                pl.col("GENOTYPE").is_in(HEMIZYGOUS_ALT_GT)
                |
                pl.col("GENOTYPE").is_in(HOM_ALT_GT)
            )
        )
        .then(6.0)

        .otherwise(0.0)

        .alias("_SCORE_INHERITANCE_COMPONENT")
    )

    # --------------------------------------------------------
    # Bounded Exomiser rank bonus.
    #
    # Exomiser rank remains important but cannot override strong
    # contradictory/weak variant evidence.
    # --------------------------------------------------------
    df = df.with_columns(
        pl.when(
            pl.col("#RANK").is_null()
        )
        .then(0.0)

        .when(
            pl.col("#RANK") <= 5
        )
        .then(12.0)

        .when(
            pl.col("#RANK") <= 10
        )
        .then(8.0)

        .when(
            pl.col("#RANK") <= 20
        )
        .then(5.0)

        .when(
            pl.col("#RANK") <= 50
        )
        .then(2.0)

        .otherwise(0.0)

        .alias("_EXOMISER_RANK_BONUS")
    )

    # --------------------------------------------------------
    # General final eligibility gate.
    #
    # A variant is eligible if it already passed KEEP_VARIANT and
    # has one of the following:
    #
    # A. strong pathogenic/curated evidence
    # B. high-impact consequence + phenotype support
    # C. coding/non-weak consequence + predictor + phenotype support
    # D. weak/noncoding consequence with specific splice or strong
    #    curated/pathogenic evidence
    #
    # This prevents generic upstream/downstream/intergenic VUSs from
    # surviving solely because their Exomiser variant score is high.
    # --------------------------------------------------------
    strong_evidence = (
        pl.col("_SCORE_CLINVAR_PLP")
        |
        pl.col("_SCORE_EXOMISER_PLP")
        |
        (pl.col("WHITELIST_VARIANT") == 1)
    )

    phenotype_supported = (
        pl.col("EXOMISER_GENE_PHENO_SCORE")
        .fill_null(0.0)
        >= 0.30
    )

    high_impact_supported = (
        pl.col("HIGH_IMPACT")
        &
        phenotype_supported
    )

    coding_predictor_supported = (
        (~pl.col("_SCORE_WEAK_NONCODING"))
        &
        pl.col("PREDICTOR_PASS")
        &
        phenotype_supported
    )

    noncoding_supported = (
        pl.col("_SCORE_WEAK_NONCODING")
        &
        (
            pl.col("_SCORE_SPLICE_SUPPORT")
            |
            strong_evidence
        )
        &
        phenotype_supported
    )

    df = df.with_columns(
        (
            pl.col("KEEP_VARIANT")
            &
            (
                strong_evidence
                |
                high_impact_supported
                |
                coding_predictor_supported
                |
                noncoding_supported
            )
        )
        .alias("FINAL_ELIGIBLE")
    )

    # --------------------------------------------------------
    # Internal score.
    #
    # Continuous evidence:
    #   phenotype score               0-28
    #   Exomiser variant score        0-16
    #   Exomiser gene combined score  0-8
    #
    # Bounded bonuses:
    #   Exomiser rank                 0-12
    #   ClinVar P/LP                  +12
    #   Exomiser P/LP                 +8
    #   high-impact                   +7
    #   predictor support             +3
    #   whitelist                     +3
    #   contributing                  +2
    #   rarity                        +0 to +4
    #   inheritance                   -6 to +7
    #
    # Penalties:
    #   ClinVar benign/conflict       -18
    #   weak/noncoding consequence    -10
    #
    # Score is clipped to 0-100.
    # --------------------------------------------------------
    df = df.with_columns(
        (
            28.0
            * pl.col("EXOMISER_GENE_PHENO_SCORE")
            .fill_null(0.0)
            .clip(0.0, 1.0)

            + 16.0
            * pl.col("EXOMISER_VARIANT_SCORE")
            .fill_null(0.0)
            .clip(0.0, 1.0)

            + 8.0
            * pl.col("EXOMISER_GENE_COMBINED_SCORE")
            .fill_null(0.0)
            .clip(0.0, 1.0)

            + pl.col("_EXOMISER_RANK_BONUS")

            + pl.when(
                pl.col("_SCORE_CLINVAR_PLP")
            ).then(12.0).otherwise(0.0)

            + pl.when(
                pl.col("_SCORE_EXOMISER_PLP")
            ).then(8.0).otherwise(0.0)

            + pl.when(
                pl.col("HIGH_IMPACT")
            ).then(7.0).otherwise(0.0)

            + pl.when(
                pl.col("PREDICTOR_PASS")
            ).then(3.0).otherwise(0.0)

            + pl.when(
                pl.col("WHITELIST_VARIANT") == 1
            ).then(3.0).otherwise(0.0)

            + pl.when(
                pl.col("CONTRIBUTING_VARIANT") == 1
            ).then(2.0).otherwise(0.0)

            + pl.col("_SCORE_FREQUENCY_BONUS")

            + pl.col("_SCORE_INHERITANCE_COMPONENT")

            - pl.when(
                pl.col("_SCORE_CLINVAR_CONFLICT")
            ).then(18.0).otherwise(0.0)

            - pl.when(
                pl.col("_SCORE_WEAK_NONCODING")
            ).then(10.0).otherwise(0.0)
        )
        .clip(0.0, 100.0)
        .round(2)
        .alias("INTERNAL_PRIORITY_SCORE")
    )

    return df


# ============================================================
# Top-N diagnostic variants
# ============================================================

def build_top_variants(df, top_n):
    return (
        df
        .filter(
            pl.col("KEEP_VARIANT")
        )
        .sort(
            by=[
                "#RANK",
                "INTERNAL_PRIORITY_SCORE",
                "WHITELIST_VARIANT",
                "CONTRIBUTING_VARIANT",
                "EXOMISER_GENE_COMBINED_SCORE",
                "EXOMISER_GENE_PHENO_SCORE",
                "EXOMISER_VARIANT_SCORE",
                "MAX_GNOMAD_SUBPOP_AF",
            ],
            descending=[
                False,
                True,
                True,
                True,
                True,
                True,
                True,
                False,
            ],
            nulls_last=True,
        )
        .unique(
            subset=VARIANT_KEY,
            keep="first",
            maintain_order=True,
        )
        .head(top_n)
    )


def format_top_output(df):
    """
    Clean and reorder the final Top-N diagnostic TSV.

    Internal analysis continues to use CONTIG/START and the prioritization
    flags. This function is applied only after ranking has finished.
    """

    # Remove duplicate/internal/debug columns without failing if a column
    # is absent in a particular run.
    df = df.drop(
        TOP_OUTPUT_DROP_COLUMNS,
        strict=False,
    )

    # Put the main variant identity and clinically useful fields first.
    first_columns = [
        column
        for column in TOP_OUTPUT_FIRST_COLUMNS
        if column in df.columns
    ]

    remaining_columns = [
        column
        for column in df.columns
        if column not in first_columns
    ]

    return df.select(
        first_columns + remaining_columns
    )

# ============================================================
# Possible compound heterozygous variants
# ============================================================

def build_compound_hets(df):

    if "EXOMISER_ACMG_DISEASE_ID" in df.columns:

        disease_id = (
            pl.col("EXOMISER_ACMG_DISEASE_ID")
            .cast(pl.String)
            .fill_null("")
        )

    else:

        disease_id = pl.lit("")

    if "EXOMISER_ACMG_DISEASE_NAME" in df.columns:

        disease_name = (
            pl.col("EXOMISER_ACMG_DISEASE_NAME")
            .cast(pl.String)
            .fill_null("")
        )

    else:

        disease_name = pl.lit("")

    df = df.with_columns(

        pl.when(
            disease_id != ""
        )
        .then(
            disease_id
        )

        .when(
            disease_name != ""
        )
        .then(
            disease_name
        )

        .otherwise(
            pl.lit("UNSPECIFIED")
        )
        .alias("_COMPOUND_DISEASE_CONTEXT")
    )

    candidates = (
        df
        .filter(
            (
                pl.col("MOI").is_in(AR_MOI)
                &
                pl.col("GENOTYPE").is_in(HET_GT)
                &
                pl.col("POPULATION_PASS")
                &
                (
                    pl.col("KEEP_VARIANT")
                    |
                    pl.col("RESCUE_VARIANT")
                    |
                    pl.col("HIGH_IMPACT")
                    |
                    pl.col("PREDICTOR_PASS")
                )
            )
        )
        .unique(
            subset=[
                "GENE_SYMBOL",
                "_COMPOUND_DISEASE_CONTEXT",
                "MOI",
                *VARIANT_KEY,
            ],
            keep="first",
            maintain_order=True,
        )
    )

    if candidates.height == 0:
        return pl.DataFrame()

    pair_rows = []

    group_columns = [
        "GENE_SYMBOL",
        "_COMPOUND_DISEASE_CONTEXT",
        "MOI",
    ]

    for group_df in candidates.partition_by(
        group_columns,
        maintain_order=True
    ):

        if group_df.height < 2:
            continue

        rows = group_df.to_dicts()

        for v1, v2 in itertools.combinations(
            rows,
            2
        ):

            key1 = tuple(
                v1.get(col)
                for col in VARIANT_KEY
            )

            key2 = tuple(
                v2.get(col)
                for col in VARIANT_KEY
            )

            if key1 == key2:
                continue

            pair_rows.append({

                "GENE_SYMBOL":
                    v1.get("GENE_SYMBOL"),

                "DISEASE_CONTEXT":
                    v1.get(
                        "_COMPOUND_DISEASE_CONTEXT"
                    ),

                "MOI":
                    v1.get("MOI"),

                "PAIR_STATUS":
                    "POSSIBLE_COMPOUND_HET",

                "PHASE_STATUS":
                    "UNKNOWN",

                "VARIANT_1_CONTIG":
                    v1.get("CONTIG"),

                "VARIANT_1_START":
                    v1.get("START"),

                "VARIANT_1_REF":
                    v1.get("REF"),

                "VARIANT_1_ALT":
                    v1.get("ALT"),

                "VARIANT_1_GT":
                    v1.get("GENOTYPE"),

                "VARIANT_1_HGVS":
                    v1.get("HGVS"),

                "VARIANT_1_FUNCTIONAL_CLASS":
                    v1.get("FUNCTIONAL_CLASS"),

                "VARIANT_1_CLINVAR":
                    v1.get(
                        "CLINVAR_PRIMARY_INTERPRETATION"
                    ),

                "VARIANT_1_EXOMISER_CLASS":
                    v1.get(
                        "EXOMISER_ACMG_CLASSIFICATION"
                    ),

                "VARIANT_1_MAX_AF":
                    v1.get(
                        "MAX_GNOMAD_SUBPOP_AF"
                    ),

                "VARIANT_2_CONTIG":
                    v2.get("CONTIG"),

                "VARIANT_2_START":
                    v2.get("START"),

                "VARIANT_2_REF":
                    v2.get("REF"),

                "VARIANT_2_ALT":
                    v2.get("ALT"),

                "VARIANT_2_GT":
                    v2.get("GENOTYPE"),

                "VARIANT_2_HGVS":
                    v2.get("HGVS"),

                "VARIANT_2_FUNCTIONAL_CLASS":
                    v2.get("FUNCTIONAL_CLASS"),

                "VARIANT_2_CLINVAR":
                    v2.get(
                        "CLINVAR_PRIMARY_INTERPRETATION"
                    ),

                "VARIANT_2_EXOMISER_CLASS":
                    v2.get(
                        "EXOMISER_ACMG_CLASSIFICATION"
                    ),

                "VARIANT_2_MAX_AF":
                    v2.get(
                        "MAX_GNOMAD_SUBPOP_AF"
                    ),
            })

    if not pair_rows:
        return pl.DataFrame()

    return pl.DataFrame(pair_rows)

# ============================================================
# Secondary findings master preparation
# ============================================================

def prepare_secondary_master(sf_df):

    required = [
        "GENE_SYMBOL",
        "INHERITANCE",
        "REPORTABLE_AS_SF",
        "VARIANTS_TO_REPORT",
    ]

    validate_columns(
        sf_df,
        required,
        "Secondary findings master file"
    )

    sf_df = sf_df.with_columns(

        pl.col("GENE_SYMBOL")
        .cast(pl.String)
        .str.strip_chars()
        .str.to_uppercase(),

        pl.col("INHERITANCE")
        .cast(pl.String)
        .str.strip_chars()
        .str.to_uppercase(),

        pl.col("REPORTABLE_AS_SF")
        .cast(pl.String)
        .str.strip_chars()
        .str.to_uppercase(),

        pl.col("VARIANTS_TO_REPORT")
        .cast(pl.String)
        .fill_null("")
        .str.strip_chars(),
    )

    sf_df = sf_df.filter(
        pl.col("REPORTABLE_AS_SF") == "YES"
    )

    if "MONDO_ID" in sf_df.columns:

        sf_df = sf_df.with_columns(

            pl.when(
                pl.col("MONDO_ID")
                .cast(pl.String)
                .fill_null("")
                != ""
            )
            .then(
                pl.col("MONDO_ID")
                .cast(pl.String)
            )
            .otherwise(
                pl.col("DISEASE_PHENOTYPE")
                .cast(pl.String)
                .fill_null("UNSPECIFIED")
                if "DISEASE_PHENOTYPE" in sf_df.columns
                else pl.lit("UNSPECIFIED")
            )
            .alias("_SF_DISEASE_KEY")
        )

    elif "DISORDER_MIM" in sf_df.columns:

        sf_df = sf_df.with_columns(

            pl.when(
                pl.col("DISORDER_MIM")
                .cast(pl.String)
                .fill_null("")
                != ""
            )
            .then(
                pl.col("DISORDER_MIM")
                .cast(pl.String)
            )
            .otherwise(
                pl.col("DISEASE_PHENOTYPE")
                .cast(pl.String)
                .fill_null("UNSPECIFIED")
                if "DISEASE_PHENOTYPE" in sf_df.columns
                else pl.lit("UNSPECIFIED")
            )
            .alias("_SF_DISEASE_KEY")
        )

    elif "DISEASE_PHENOTYPE" in sf_df.columns:

        sf_df = sf_df.with_columns(

            pl.col("DISEASE_PHENOTYPE")
            .cast(pl.String)
            .fill_null("UNSPECIFIED")
            .alias("_SF_DISEASE_KEY")
        )

    else:

        sf_df = sf_df.with_columns(
            pl.lit("UNSPECIFIED")
            .alias("_SF_DISEASE_KEY")
        )

    return sf_df

# ============================================================
# Secondary findings
# ============================================================

def build_secondary_findings(
    variant_df,
    secondary_gene_df
):

    sf_master = prepare_secondary_master(
        secondary_gene_df
    )

    print()
    print("========== SECONDARY FINDINGS DEBUG ==========")

    print(
        f"[SF] SF master rows              : "
        f"{sf_master.height}"
    )

    print(
        f"[SF] Unique SF genes             : "
        f"{sf_master['GENE_SYMBOL'].n_unique()}"
    )

    input_genes = set(
        variant_df
        .get_column("GENE_SYMBOL")
        .drop_nulls()
        .unique()
        .to_list()
    )

    sf_genes = set(
        sf_master
        .get_column("GENE_SYMBOL")
        .drop_nulls()
        .unique()
        .to_list()
    )

    overlap = input_genes & sf_genes

    print(
        f"[SF] Input unique genes          : "
        f"{len(input_genes)}"
    )

    print(
        f"[SF] SF genes present in input   : "
        f"{len(overlap)}"
    )

    if overlap:

        print(
            "[SF] Example overlapping genes   : "
            + ", ".join(
                sorted(overlap)[:20]
            )
        )

    # --------------------------------------------------------
    # Gene match only
    #
    # No HPO / Exomiser rank / CONTRIBUTING filter.
    # --------------------------------------------------------

    secondary = variant_df.join(
        sf_master,
        on="GENE_SYMBOL",
        how="inner",
        suffix="_SF",
    )

    print(
        f"[SF] Rows after gene join        : "
        f"{secondary.height}"
    )

    if secondary.height == 0:

        print(
            "[SF] STOP: no variants overlap "
            "with SF genes."
        )

        return pl.DataFrame()

    # --------------------------------------------------------
    # P/LP candidate evidence
    # --------------------------------------------------------

    secondary = secondary.with_columns(

        plp_expr(
            "CLINVAR_PRIMARY_INTERPRETATION"
        ).alias(
            "_CLINVAR_PLP"
        ),

        plp_expr(
            "EXOMISER_ACMG_CLASSIFICATION"
        ).alias(
            "_EXOMISER_PLP"
        ),
    )

    clinvar_count = secondary.filter(
        pl.col("_CLINVAR_PLP")
    ).height

    exomiser_count = secondary.filter(
        pl.col("_EXOMISER_PLP")
    ).height

    print(
        f"[SF] ClinVar P/LP rows           : "
        f"{clinvar_count}"
    )

    print(
        f"[SF] Exomiser P/LP rows          : "
        f"{exomiser_count}"
    )

    secondary = secondary.with_columns(

        (
            pl.col("_CLINVAR_PLP")
            |
            pl.col("_EXOMISER_PLP")
        )
        .alias("_P_LP_CANDIDATE")
    )

    secondary = secondary.filter(
        pl.col("_P_LP_CANDIDATE")
    )

    print(
        f"[SF] Rows after P/LP filter      : "
        f"{secondary.height}"
    )

    if secondary.height == 0:

        print(
            "[SF] STOP: SF genes exist, but "
            "no P/LP candidates were found."
        )

        return pl.DataFrame()

    # --------------------------------------------------------
    # Source
    # --------------------------------------------------------

    secondary = secondary.with_columns(

        pl.when(
            pl.col("_CLINVAR_PLP")
            &
            pl.col("_EXOMISER_PLP")
        )
        .then(
            pl.lit("CLINVAR+EXOMISER")
        )

        .when(
            pl.col("_CLINVAR_PLP")
        )
        .then(
            pl.lit("CLINVAR")
        )

        .when(
            pl.col("_EXOMISER_PLP")
        )
        .then(
            pl.lit("EXOMISER")
        )

        .otherwise(
            pl.lit(None)
        )
        .alias("SF_CLASSIFICATION_SOURCE")
    )

    # --------------------------------------------------------
    # Genotype
    # --------------------------------------------------------

    secondary = secondary.with_columns(

        pl.col("GENOTYPE")
        .is_in(HET_GT)
        .alias("_IS_HET"),

        pl.col("GENOTYPE")
        .is_in(HOM_ALT_GT)
        .alias("_IS_HOM_ALT"),

        pl.col("GENOTYPE")
        .is_in(HEMIZYGOUS_ALT_GT)
        .alias("_IS_HEMIZYGOUS"),
    )

    sf_group = [
        "GENE_SYMBOL",
        "_SF_DISEASE_KEY",
        "INHERITANCE",
    ]

    # --------------------------------------------------------
    # Count distinct heterozygous variants per SF model
    # --------------------------------------------------------

    secondary = secondary.with_columns(

        pl.when(
            pl.col("_IS_HET")
        )
        .then(
            pl.concat_str(
                [
                    pl.col("CONTIG"),
                    pl.col("START").cast(pl.String),
                    pl.col("REF"),
                    pl.col("ALT"),
                ],
                separator=":"
            )
        )
        .otherwise(
            pl.lit(None)
        )
        .alias("_HET_VARIANT_KEY")
    )

    het_counts = (
        secondary
        .filter(
            pl.col("_IS_HET")
        )
        .group_by(
            sf_group
        )
        .agg(
            pl.col("_HET_VARIANT_KEY")
            .n_unique()
            .alias("SF_HET_PLP_COUNT")
        )
    )

    secondary = secondary.join(
        het_counts,
        on=sf_group,
        how="left",
    )

    secondary = secondary.with_columns(

        pl.col("SF_HET_PLP_COUNT")
        .fill_null(0)
    )

    # --------------------------------------------------------
    # Detect wording requiring 2 het variants
    # --------------------------------------------------------

    secondary = secondary.with_columns(

        pl.col("VARIANTS_TO_REPORT")
        .str.to_uppercase()
        .str.contains(
            r"2\s*HET|"
            r"TWO\s+HET|"
            r"2\s*HETEROZYGOUS|"
            r"TWO\s+HETEROZYGOUS"
        )
        .fill_null(False)
        .alias("_REQUIRES_TWO_HET")
    )

    # --------------------------------------------------------
    # Allelic state
    # --------------------------------------------------------

    secondary = secondary.with_columns(

        pl.when(
            pl.col("_IS_HEMIZYGOUS")
        )
        .then(
            pl.lit("HEMIZYGOUS")
        )

        .when(
            pl.col("_IS_HOM_ALT")
        )
        .then(
            pl.lit("HOMOZYGOUS_ALT")
        )

        .when(
            pl.col("_IS_HET")
        )
        .then(
            pl.lit("HETEROZYGOUS")
        )

        .otherwise(
            pl.lit("OTHER")
        )
        .alias("SF_ALLELIC_STATE")
    )

    # --------------------------------------------------------
    # SF candidate status
    # --------------------------------------------------------

    secondary = secondary.with_columns(

        # AD
        pl.when(
            pl.col("INHERITANCE").is_in(AD_MOI)
            &
            (
                pl.col("_IS_HET")
                |
                pl.col("_IS_HOM_ALT")
            )
        )
        .then(
            pl.lit("SF_CANDIDATE_REVIEW")
        )

        # AR homozygous
        .when(
            pl.col("INHERITANCE").is_in(AR_MOI)
            &
            pl.col("_IS_HOM_ALT")
        )
        .then(
            pl.lit("SF_CANDIDATE_REVIEW")
        )

        # AR possible compound het
        .when(
            pl.col("INHERITANCE").is_in(AR_MOI)
            &
            pl.col("_IS_HET")
            &
            (
                pl.col("SF_HET_PLP_COUNT") >= 2
            )
        )
        .then(
            pl.lit("PHASING_REQUIRED")
        )

        # XL hemi/hom
        .when(
            pl.col("INHERITANCE").is_in(XL_MOI)
            &
            (
                pl.col("_IS_HEMIZYGOUS")
                |
                pl.col("_IS_HOM_ALT")
            )
        )
        .then(
            pl.lit("SF_CANDIDATE_REVIEW")
        )

        # XL requiring two hets
        .when(
            pl.col("INHERITANCE").is_in(XL_MOI)
            &
            pl.col("_IS_HET")
            &
            pl.col("_REQUIRES_TWO_HET")
            &
            (
                pl.col("SF_HET_PLP_COUNT") >= 2
            )
        )
        .then(
            pl.lit("PHASING_REQUIRED")
        )

        # XL heterozygous and rule does not explicitly require 2 het
        .when(
            pl.col("INHERITANCE").is_in(XL_MOI)
            &
            pl.col("_IS_HET")
            &
            (~pl.col("_REQUIRES_TWO_HET"))
        )
        .then(
            pl.lit("SF_CANDIDATE_REVIEW")
        )

        # Unknown / uncommon MOI
        .when(
            ~pl.col("INHERITANCE").is_in(
                AD_MOI + AR_MOI + XL_MOI
            )
        )
        .then(
            pl.lit("SF_MANUAL_REVIEW")
        )

        .otherwise(
            pl.lit("NOT_REPORTABLE_ALLELIC_STATE")
        )
        .alias("SF_STATUS")
    )

    secondary = secondary.with_columns(

        pl.when(
            pl.col("SF_STATUS") == "PHASING_REQUIRED"
        )
        .then(
            pl.lit("POSSIBLE_BIALLELIC")
        )
        .otherwise(
            pl.col("SF_ALLELIC_STATE")
        )
        .alias("SF_ALLELIC_STATE"),

        pl.when(
            pl.col("SF_STATUS") == "PHASING_REQUIRED"
        )
        .then(
            pl.lit("UNKNOWN")
        )
        .otherwise(
            pl.lit(None)
        )
        .alias("SF_PHASE_STATUS"),
    )

    print(
        "[SF] Status counts:"
    )

    print(
        secondary
        .group_by("SF_STATUS")
        .len()
        .sort("len", descending=True)
    )

    # --------------------------------------------------------
    # Keep candidates/reviewable findings
    # --------------------------------------------------------

    secondary = secondary.filter(

        pl.col("SF_STATUS").is_in(
            [
                "SF_CANDIDATE_REVIEW",
                "PHASING_REQUIRED",
                "SF_MANUAL_REVIEW",
            ]
        )
    )

    print(
        f"[SF] Final SF candidate rows     : "
        f"{secondary.height}"
    )

    print(
        "=============================================="
    )
    print()

    if secondary.height == 0:
        return pl.DataFrame()

    # --------------------------------------------------------
    # Deduplicate disease-specific physical variants
    # --------------------------------------------------------

    secondary = secondary.unique(

        subset=[
            "GENE_SYMBOL",
            "_SF_DISEASE_KEY",
            *VARIANT_KEY,
        ],

        keep="first",
        maintain_order=True,
    )

    secondary = secondary.sort(

        by=[
            "GENE_SYMBOL",
            "_SF_DISEASE_KEY",
            "SF_STATUS",
            "CONTIG",
            "START",
        ],

        nulls_last=True,
    )

    secondary = secondary.drop(
        [
            "_CLINVAR_PLP",
            "_EXOMISER_PLP",
            "_P_LP_CANDIDATE",
            "_IS_HET",
            "_IS_HOM_ALT",
            "_IS_HEMIZYGOUS",
            "_REQUIRES_TWO_HET",
            "_HET_VARIANT_KEY",
            "_SF_DISEASE_KEY",
        ],
        strict=False,
    )

    return secondary

# ============================================================
# Output writer
# ============================================================

def write_dataframe(df, path):

    if df.width == 0:

        path.write_text(
            "",
            encoding="utf-8"
        )

    else:

        df.write_csv(
            path,
            separator="\t"
        )

# ============================================================
# Main
# ============================================================

def main():

    args = parse_arguments()

    input_path = Path(
        args.input
    )

    secondary_gene_path = Path(
        args.secondary_gene_list
    )

    secondary_input_path = (
        Path(args.secondary_input)
        if args.secondary_input
        else input_path
    )

    output_dir = Path(
        args.output_dir
    )

    try:

        if not input_path.is_file():

            raise FileNotFoundError(
                f"Input TSV not found: "
                f"{input_path}"
            )

        if not secondary_gene_path.is_file():

            raise FileNotFoundError(
                f"Secondary findings master not found: "
                f"{secondary_gene_path}"
            )

        if not secondary_input_path.is_file():

            raise FileNotFoundError(
                f"Secondary input TSV not found: "
                f"{secondary_input_path}"
            )

        output_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        # ====================================================
        # Diagnostic input
        # ====================================================

        df = read_tsv(
            input_path
        )

        required_columns = [
            "#RANK",
            "GENE_SYMBOL",
            "MOI",
            "CONTIG",
            "START",
            "REF",
            "ALT",
            "GENOTYPE",
            "FUNCTIONAL_CLASS",
            "CONTRIBUTING_VARIANT",
            "WHITELIST_VARIANT",
            "EXOMISER_GENE_COMBINED_SCORE",
            "EXOMISER_GENE_PHENO_SCORE",
            "EXOMISER_VARIANT_SCORE",
            "EXOMISER_ACMG_CLASSIFICATION",
            "CLINVAR_PRIMARY_INTERPRETATION",
            "CADD_RAW_SCORE",
            "REVEL_SCORE",
            "SPLICE_AI_SCORE",
        ]

        validate_columns(
            df,
            required_columns,
            "Merged Exomiser TSV"
        )

        df = prepare_dataframe(
            df
        )

        df = add_prioritization_flags(
            df
        )

        # ====================================================
        # Internal priority score
        # ====================================================

        df = add_internal_priority_score(
            df
        )

        # ====================================================
        # Top-N
        # ====================================================

        top_df = build_top_variants(
            df,
            args.top
        )

        # Clean/reorder only the final Top-N output.
        # The full internal dataframe remains unchanged for compound-het
        # and secondary-finding analysis.
        top_df = format_top_output(
            top_df
        )

        # ====================================================
        # Compound het
        # ====================================================

        compound_df = build_compound_hets(
            df
        )

        # ====================================================
        # Secondary findings master
        # ====================================================

        secondary_gene_df = read_tsv(
            secondary_gene_path
        )

        # ====================================================
        # Secondary findings variant source
        #
        # Prefer full QC-passed annotated WGS table.
        # ====================================================

        if secondary_input_path == input_path:

            secondary_variant_df = df.clone()

            print(
                "WARNING: --secondary-input was not supplied."
            )

            print(
                "Secondary findings will be searched only in "
                "the main Exomiser-merged input."
            )

        else:

            secondary_variant_df = read_tsv(
                secondary_input_path
            )

            secondary_variant_df = prepare_secondary_input(
                secondary_variant_df
            )

        # ====================================================
        # Secondary findings
        # ====================================================

        secondary_df = build_secondary_findings(
            secondary_variant_df,
            secondary_gene_df
        )

        # ====================================================
        # Outputs
        # ====================================================

        top_output = (
            output_dir
            / f"{args.prefix}_top{args.top}.tsv"
        )

        compound_output = (
            output_dir
            / f"{args.prefix}_possible_compound_het.tsv"
        )

        secondary_output = (
            output_dir
            / f"{args.prefix}_secondary_findings.tsv"
        )

        write_dataframe(
            top_df,
            top_output
        )

        write_dataframe(
            compound_df,
            compound_output
        )

        write_dataframe(
            secondary_df,
            secondary_output
        )

        # ====================================================
        # Summary
        # ====================================================

        print()
        print("============== SUMMARY ==============")

        print(
            f"Diagnostic input rows          : "
            f"{df.height}"
        )

        print(
            f"Prioritized Top-{args.top} rows      : "
            f"{top_df.height}"
        )

        print(
            f"Possible compound-het pairs    : "
            f"{compound_df.height}"
        )

        print(
            f"Secondary source rows          : "
            f"{secondary_variant_df.height}"
        )

        print(
            f"Secondary finding candidates   : "
            f"{secondary_df.height}"
        )

        print()

        print(
            f"Top variants                   : "
            f"{top_output}"
        )

        print(
            f"Compound hets                  : "
            f"{compound_output}"
        )

        print(
            f"Secondary findings             : "
            f"{secondary_output}"
        )

        print("=====================================")

        return 0

    except Exception as error:

        print(
            f"ERROR: {error}",
            file=sys.stderr
        )

        return 1

if __name__ == "__main__":
    sys.exit(
        main()
    )