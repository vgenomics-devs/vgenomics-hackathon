#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import polars as pl

JOIN_KEYS = [
    "CONTIG",
    "START",
    "REF",
    "ALT",
]

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Merge Exomiser variants TSV with Exomiser JSON-derived TSV "
            "using CONTIG, START, REF and ALT while preserving the "
            "original variants TSV row order."
        )
    )
    parser.add_argument(
        "--variants-tsv",
        required=True,
        help="Exomiser variants TSV file."
    )
    parser.add_argument(
        "--json-tsv",
        required=True,
        help="TSV extracted from Exomiser JSON output."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output merged TSV file."
    )
    return parser.parse_args()

def validate_columns(
    df: pl.DataFrame,
    required_columns: list[str],
    label: str,
):
    missing = [
        column
        for column in required_columns
        if column not in df.columns
    ]
    if missing:
        raise ValueError(
            f"{label} is missing required column(s): "
            + ", ".join(missing)
        )

def normalize_join_columns(
    df: pl.DataFrame,
) -> pl.DataFrame:
    return df.with_columns(
        pl.col("CONTIG")
        .cast(pl.String)
        .str.strip_chars()
        .str.replace(r"^chr", "")
        .replace({
            "M": "MT",
        })
        .alias("CONTIG"),
        pl.col("START")
        .cast(pl.Int64, strict=False)
        .alias("START"),
        pl.col("REF")
        .cast(pl.String)
        .str.strip_chars()
        .str.to_uppercase()
        .alias("REF"),
        pl.col("ALT")
        .cast(pl.String)
        .str.strip_chars()
        .str.to_uppercase()
        .alias("ALT"),
    )

def main():
    args = parse_arguments()
    variants_path = Path(args.variants_tsv)
    json_path = Path(args.json_tsv)
    output_path = Path(args.output)
    try:
        # ====================================================
        # Validate input files
        # ====================================================
        if not variants_path.is_file():
            raise FileNotFoundError(
                f"Variants TSV not found: {variants_path}"
            )
        if not json_path.is_file():
            raise FileNotFoundError(
                f"JSON TSV not found: {json_path}"
            )
        output_path.parent.mkdir(
            parents=True,
            exist_ok=True
        )
        # ====================================================
        # Read TSV files
        # ====================================================
        variants_df = pl.read_csv(
            variants_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[".", ""],
        )
        json_df = pl.read_csv(
            json_path,
            separator="\t",
            infer_schema_length=10000,
            null_values=[".", ""],
        )
        # ====================================================
        # Validate join columns
        # ====================================================
        validate_columns(
            variants_df,
            JOIN_KEYS,
            "Exomiser variants TSV",
        )
        validate_columns(
            json_df,
            JOIN_KEYS,
            "Exomiser JSON TSV",
        )
        # ====================================================
        # Add explicit row-order column
        #
        # This makes order preservation deterministic even if
        # the join implementation changes in future.
        # ====================================================
        variants_df = variants_df.with_row_index(
            "_VARIANT_ORDER"
        )
        # ====================================================
        # Normalize join keys
        # ====================================================
        variants_df = normalize_join_columns(
            variants_df
        )
        json_df = normalize_join_columns(
            json_df
        )
        # ====================================================
        # Validate JSON TSV uniqueness
        #
        # JSON TSV should contain one row per physical variant.
        # If duplicate keys exist, a join would duplicate rows
        # from variants.tsv.
        # ====================================================
        duplicate_json_keys = (
            json_df
            .group_by(JOIN_KEYS)
            .len()
            .filter(
                pl.col("len") > 1
            )
        )
        if duplicate_json_keys.height > 0:
            raise ValueError(
                f"JSON TSV contains "
                f"{duplicate_json_keys.height} duplicate "
                "CONTIG/START/REF/ALT variant key(s). "
                "Expected one row per physical variant."
            )
        # ====================================================
        # Remove redundant columns from JSON TSV
        #
        # Join keys already exist in variants.tsv.
        # Only annotation columns from JSON TSV are appended.
        # ====================================================
        json_annotation_columns = [
            column
            for column in json_df.columns
            if column not in JOIN_KEYS
        ]
        # ====================================================
        # Avoid duplicate non-key column names
        #
        # If a JSON column already exists in variants.tsv,
        # rename the JSON version with JSON_ prefix.
        # ====================================================
        rename_map = {}
        for column in json_annotation_columns:
            if column in variants_df.columns:
                rename_map[column] = (
                    f"JSON_{column}"
                )
        if rename_map:
            json_df = json_df.rename(
                rename_map
            )
        # ====================================================
        # Merge
        #
        # LEFT JOIN:
        # - every variants.tsv row is retained
        # - JSON annotations are appended where matched
        # ====================================================
        merged_df = variants_df.join(
            json_df,
            on=JOIN_KEYS,
            how="left",
        )
        # ====================================================
        # Restore exact variants.tsv order
        # ====================================================
        merged_df = (
            merged_df
            .sort("_VARIANT_ORDER")
            .drop("_VARIANT_ORDER")
        )
        # ====================================================
        # Write output
        # ====================================================
        merged_df.write_csv(
            output_path,
            separator="\t",
        )
        # ====================================================
        # Summary
        # ====================================================
        print(
            f"Variants TSV rows : {variants_df.height}"
        )
        print(
            f"JSON TSV rows     : {json_df.height}"
        )
        print(
            f"Merged rows       : {merged_df.height}"
        )
        print(
            f"Output            : {output_path}"
        )
        return 0
    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        return 1

if __name__ == "__main__":
    sys.exit(
        main()
    )