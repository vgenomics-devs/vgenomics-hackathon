#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import polars as pl


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Merge VCF-extracted TSV with enriched Exomiser TSV using "
            "VCF CHROM/POS/REF/ALT against Exomiser CONTIG/START/REF/ALT. "
            "Final output uses canonical CONTIG/START/END/REF/ALT columns."
        )
    )
    parser.add_argument("--vcf-tsv", required=True)
    parser.add_argument("--exomiser-tsv", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def normalize_chromosome_expr(column_name):
    return (
        pl.col(column_name)
        .cast(pl.String)
        .str.strip_chars()
        .str.replace(r"^chr", "")
        .replace({"M": "MT"})
    )


def normalize_allele_expr(column_name):
    return (
        pl.col(column_name)
        .cast(pl.String)
        .str.strip_chars()
        .str.to_uppercase()
    )


def validate_columns(df, required_columns, file_label):
    missing_columns = [
        column for column in required_columns
        if column not in df.columns
    ]
    if missing_columns:
        raise ValueError(
            f"{file_label} is missing required column(s): "
            + ", ".join(missing_columns)
        )


def main():
    args = parse_arguments()

    vcf_path = Path(args.vcf_tsv)
    exomiser_path = Path(args.exomiser_tsv)
    output_path = Path(args.output)

    try:
        if not vcf_path.is_file():
            raise FileNotFoundError(f"VCF TSV not found: {vcf_path}")

        if not exomiser_path.is_file():
            raise FileNotFoundError(f"Exomiser TSV not found: {exomiser_path}")

        output_path.parent.mkdir(parents=True, exist_ok=True)

        vcf_df = pl.read_csv(
            vcf_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[".", ""],
            schema_overrides={
                "CHROM": pl.String,
                "POS": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
                "GT": pl.String,
                "ZYGOSITY": pl.String,
                "AD": pl.String,
            }
        )

        exomiser_df = pl.read_csv(
            exomiser_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[".", ""],
            schema_overrides={
                "CONTIG": pl.String,
                "START": pl.Int64,
                "END": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
            }
        )

        validate_columns(
            vcf_df,
            ["CHROM", "POS", "REF", "ALT"],
            "VCF TSV"
        )

        validate_columns(
            exomiser_df,
            ["CONTIG", "START", "REF", "ALT"],
            "Exomiser TSV"
        )

        vcf_df = vcf_df.with_columns(
            normalize_chromosome_expr("CHROM").alias("_JOIN_CHROM"),
            pl.col("POS").cast(pl.Int64, strict=False).alias("_JOIN_POS"),
            normalize_allele_expr("REF").alias("_JOIN_REF"),
            normalize_allele_expr("ALT").alias("_JOIN_ALT"),
        )

        exomiser_df = exomiser_df.with_columns(
            normalize_chromosome_expr("CONTIG").alias("_JOIN_CHROM"),
            pl.col("START").cast(pl.Int64, strict=False).alias("_JOIN_POS"),
            normalize_allele_expr("REF").alias("_JOIN_REF"),
            normalize_allele_expr("ALT").alias("_JOIN_ALT"),
            pl.lit("MATCHED").alias("_EXOMISER_MATCH_MARKER"),
        )

        null_vcf_keys = (
            vcf_df
            .filter(
                pl.any_horizontal([
                    pl.col("_JOIN_CHROM").is_null(),
                    pl.col("_JOIN_POS").is_null(),
                    pl.col("_JOIN_REF").is_null(),
                    pl.col("_JOIN_ALT").is_null(),
                ])
            )
            .height
        )
        if null_vcf_keys > 0:
            raise ValueError(
                f"VCF TSV contains {null_vcf_keys} row(s) with incomplete merge keys."
            )

        null_exomiser_keys = (
            exomiser_df
            .filter(
                pl.any_horizontal([
                    pl.col("_JOIN_CHROM").is_null(),
                    pl.col("_JOIN_POS").is_null(),
                    pl.col("_JOIN_REF").is_null(),
                    pl.col("_JOIN_ALT").is_null(),
                ])
            )
            .height
        )
        if null_exomiser_keys > 0:
            raise ValueError(
                f"Exomiser TSV contains {null_exomiser_keys} row(s) with incomplete merge keys."
            )

        duplicate_vcf_keys = (
            vcf_df
            .group_by([
                "_JOIN_CHROM",
                "_JOIN_POS",
                "_JOIN_REF",
                "_JOIN_ALT"
            ])
            .len()
            .filter(pl.col("len") > 1)
        )

        if duplicate_vcf_keys.height > 0:
            raise ValueError(
                f"VCF TSV contains {duplicate_vcf_keys.height} duplicate physical variant key(s)."
            )

        duplicate_exomiser_keys = (
            exomiser_df
            .group_by([
                "_JOIN_CHROM",
                "_JOIN_POS",
                "_JOIN_REF",
                "_JOIN_ALT"
            ])
            .len()
            .filter(pl.col("len") > 1)
        )

        if duplicate_exomiser_keys.height > 0:
            print(
                f"WARNING: Exomiser TSV contains "
                f"{duplicate_exomiser_keys.height} duplicate physical variant key(s). "
                "Keeping first row for each key.",
                file=sys.stderr
            )
            exomiser_df = exomiser_df.unique(
                subset=[
                    "_JOIN_CHROM",
                    "_JOIN_POS",
                    "_JOIN_REF",
                    "_JOIN_ALT"
                ],
                keep="first",
                maintain_order=True
            )

        exomiser_drop_columns = {"CONTIG", "START", "REF", "ALT"}

        exomiser_annotation_columns = [
            column
            for column in exomiser_df.columns
            if column not in exomiser_drop_columns
        ]

        exomiser_df = exomiser_df.select(exomiser_annotation_columns)

        internal_columns = {
            "_JOIN_CHROM",
            "_JOIN_POS",
            "_JOIN_REF",
            "_JOIN_ALT",
            "_EXOMISER_MATCH_MARKER"
        }

        rename_map = {}
        for column in exomiser_df.columns:
            if column in vcf_df.columns and column not in internal_columns:
                rename_map[column] = f"EXOMISER_{column}"

        if rename_map:
            exomiser_df = exomiser_df.rename(rename_map)

        merged_df = vcf_df.join(
            exomiser_df,
            on=[
                "_JOIN_CHROM",
                "_JOIN_POS",
                "_JOIN_REF",
                "_JOIN_ALT"
            ],
            how="left"
        )

        merged_df = merged_df.with_columns(
            pl.when(pl.col("_EXOMISER_MATCH_MARKER").is_not_null())
            .then(pl.lit("MATCHED"))
            .otherwise(pl.lit("NOT_MATCHED"))
            .alias("EXOMISER_MATCH_STATUS")
        )

        merged_df = merged_df.with_columns(
            normalize_chromosome_expr("CHROM").alias("CONTIG"),
            pl.col("POS").cast(pl.Int64, strict=False).alias("START"),
            normalize_allele_expr("REF").alias("REF"),
            normalize_allele_expr("ALT").alias("ALT"),
        )

        if "END" in merged_df.columns:
            merged_df = merged_df.with_columns(
                pl.coalesce([
                    pl.col("END").cast(pl.Int64, strict=False),
                    (
                        pl.col("START")
                        + pl.col("REF").str.len_chars()
                        - 1
                    )
                ]).alias("END")
            )
        else:
            merged_df = merged_df.with_columns(
                (
                    pl.col("START")
                    + pl.col("REF").str.len_chars()
                    - 1
                ).alias("END")
            )

        merged_df = merged_df.drop([
            "_JOIN_CHROM",
            "_JOIN_POS",
            "_JOIN_REF",
            "_JOIN_ALT",
            "_EXOMISER_MATCH_MARKER"
        ])

        preferred_columns = [
            "CONTIG",
            "START",
            "END",
            "REF",
            "ALT",
            "CHROM",
            "POS",
            "QUAL",
            "FILTER",
            "INFO",
            "GT",
            "ZYGOSITY",
            "AD",
            "VAF",
            "DP",
            "GQ",
            "PL",
            "EXOMISER_MATCH_STATUS"
        ]

        existing_preferred_columns = [
            column
            for column in preferred_columns
            if column in merged_df.columns
        ]

        remaining_columns = [
            column
            for column in merged_df.columns
            if column not in existing_preferred_columns
        ]

        merged_df = merged_df.select(
            existing_preferred_columns + remaining_columns
        )

        merged_df.write_csv(output_path, separator="\t")

        total_vcf = vcf_df.height
        total_exomiser = exomiser_df.height
        matched = (
            merged_df
            .filter(pl.col("EXOMISER_MATCH_STATUS") == "MATCHED")
            .height
        )
        unmatched = total_vcf - matched

        print(f"VCF variants        : {total_vcf}")
        print(f"Exomiser variants   : {total_exomiser}")
        print(f"Matched variants    : {matched}")
        print(f"VCF-only variants   : {unmatched}")
        print(f"Output rows         : {merged_df.height}")
        print(f"Output file         : {output_path}")

        if matched == 0:
            print(
                "WARNING: No VCF variants matched Exomiser. "
                "Check chromosome naming, coordinates, REF/ALT representation "
                "and normalization.",
                file=sys.stderr
            )

            print("VCF join-key examples:", file=sys.stderr)
            print(
                vcf_df.select([
                    "_JOIN_CHROM",
                    "_JOIN_POS",
                    "_JOIN_REF",
                    "_JOIN_ALT"
                ]).head(5),
                file=sys.stderr
            )

            print("Exomiser join-key examples:", file=sys.stderr)
            print(
                exomiser_df.select([
                    "_JOIN_CHROM",
                    "_JOIN_POS",
                    "_JOIN_REF",
                    "_JOIN_ALT"
                ]).head(5),
                file=sys.stderr
            )

        return 0

    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
