#!/usr/bin/env python3
import argparse
import json
import re
import sys
from pathlib import Path
import polars as pl
# ============================================================
# Arguments
# ============================================================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Extract selected Exomiser variant-level frequency "
            "and pathogenicity information from JSON-lines output."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input Exomiser JSON/JSONL file."
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file."
    )
    return parser.parse_args()
# ============================================================
# Normalize source names
# ============================================================
def normalize_name(value):
    if value is None:
        return None
    value = str(value).strip().upper()
    value = re.sub(
        r"[^A-Z0-9]+",
        "_",
        value
    )
    return value.strip("_")
# ============================================================
# Read JSON-lines
# ============================================================
def read_exomiser_json(path):
    records = []
    with open(
        path,
        "r",
        encoding="utf-8"
    ) as handle:
        for line_number, line in enumerate(
            handle,
            start=1
        ):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Invalid JSON at line "
                    f"{line_number}: {error}"
                )
            records.append(record)
    if not records:
        raise ValueError(
            "No JSON records found in input file."
        )
    return records
# ============================================================
# Convert lists to comma-separated strings
# ============================================================
def list_to_string(values):
    if not values:
        return None
    return ",".join(
        str(value)
        for value in values
    )
# ============================================================
# Set value only if not already present
# ============================================================
def set_value(row, column, value):
    if value is None:
        return
    if column not in row or row[column] is None:
        row[column] = value
# ============================================================
# Extract variant details
# ============================================================
def extract_variants(gene_records):
    variants = {}
    frequency_sources = set()
    pathogenicity_sources = set()
    # --------------------------------------------------------
    # First pass:
    # Discover available frequency and pathogenicity sources
    # --------------------------------------------------------
    for gene_record in gene_records:
        for variant in gene_record.get(
            "variantEvaluations",
            []
        ):
            # Frequency sources
            frequency_data = (
                variant.get("frequencyData")
                or {}
            )
            for frequency in (
                frequency_data.get(
                    "frequencies"
                )
                or []
            ):
                source = normalize_name(
                    frequency.get("source")
                )
                if source:
                    frequency_sources.add(
                        source
                    )
            # Pathogenicity sources
            pathogenicity_data = (
                variant.get(
                    "pathogenicityData"
                )
                or {}
            )
            for prediction in (
                pathogenicity_data.get(
                    "pathogenicityScores"
                )
                or []
            ):
                source = normalize_name(
                    prediction.get("source")
                )
                if source:
                    pathogenicity_sources.add(
                        source
                    )
    # --------------------------------------------------------
    # Second pass:
    # One row per physical variant
    # --------------------------------------------------------
    for gene_record in gene_records:
        gene_symbol = gene_record.get(
            "geneSymbol"
        )
        for variant in gene_record.get(
            "variantEvaluations",
            []
        ):
            contig = variant.get(
                "contigName"
            )
            start = variant.get(
                "start"
            )
            ref = variant.get(
                "ref"
            )
            alt = variant.get(
                "alt"
            )
            if (
                contig is None
                or start is None
                or ref is None
                or alt is None
            ):
                continue
            key = (
                str(contig),
                int(start),
                str(ref).upper(),
                str(alt).upper(),
            )
            # ------------------------------------------------
            # Initialize row
            # ------------------------------------------------
            if key not in variants:
                variants[key] = {
                    "CONTIG":
                        str(contig),
                    "START":
                        int(start),
                    "REF":
                        str(ref).upper(),
                    "ALT":
                        str(alt).upper(),
                    "EXOMISER_JSON_GENES":
                        set(),
                    "FREQUENCY_SCORE":
                        variant.get(
                            "frequencyScore"
                        ),
                    "PATHOGENICITY_SCORE":
                        variant.get(
                            "pathogenicityScore"
                        ),
                    "PREDICTED_PATHOGENIC":
                        variant.get(
                            "predictedPathogenic"
                        ),
                    "PASSED_FILTER_TYPES":
                        list_to_string(
                            variant.get(
                                "passedFilterTypes"
                            )
                        ),
                    "COMPATIBLE_INHERITANCE_MODES":
                        list_to_string(
                            variant.get(
                                "compatibleInheritanceModes"
                            )
                        ),
                    "CONTRIBUTING_INHERITANCE_MODES":
                        list_to_string(
                            variant.get(
                                "contributingInheritanceModes"
                            )
                        ),
                }
            row = variants[key]
            if gene_symbol:
                row[
                    "EXOMISER_JSON_GENES"
                ].add(
                    gene_symbol
                )
            # =================================================
            # Frequency data
            #
            # Keep AF only.
            #
            # Do NOT extract:
            # AC
            # AN
            # HOMS
            # =================================================
            frequency_data = (
                variant.get("frequencyData")
                or {}
            )
            for frequency in (
                frequency_data.get(
                    "frequencies"
                )
                or []
            ):
                source = normalize_name(
                    frequency.get("source")
                )
                if not source:
                    continue
                set_value(
                    row,
                    f"{source}_AF",
                    frequency.get(
                        "frequency"
                    )
                )
            # =================================================
            # Pathogenicity scores
            # =================================================
            pathogenicity_data = (
                variant.get(
                    "pathogenicityData"
                )
                or {}
            )
            for prediction in (
                pathogenicity_data.get(
                    "pathogenicityScores"
                )
                or []
            ):
                source = normalize_name(
                    prediction.get("source")
                )
                if not source:
                    continue
                # --------------------------------------------
                # Generic score
                # --------------------------------------------
                set_value(
                    row,
                    f"{source}_SCORE",
                    prediction.get(
                        "score"
                    )
                )
                # --------------------------------------------
                # Raw score
                # --------------------------------------------
                set_value(
                    row,
                    f"{source}_RAW_SCORE",
                    prediction.get(
                        "rawScore"
                    )
                )
                # --------------------------------------------
                # Scaled score
                # --------------------------------------------
                set_value(
                    row,
                    f"{source}_SCALED_SCORE",
                    prediction.get(
                        "scaledScore"
                    )
                )
    # --------------------------------------------------------
    # Convert gene sets to strings
    # --------------------------------------------------------
    output_rows = []
    for row in variants.values():
        genes = row.get(
            "EXOMISER_JSON_GENES"
        )
        if isinstance(
            genes,
            set
        ):
            row[
                "EXOMISER_JSON_GENES"
            ] = ",".join(
                sorted(genes)
            )
        output_rows.append(
            row
        )
    return (
        output_rows,
        sorted(
            frequency_sources
        ),
        sorted(
            pathogenicity_sources
        ),
    )
# ============================================================
# Create Polars DataFrame
# ============================================================
def create_dataframe(
    rows,
    frequency_sources,
    pathogenicity_sources
):
    if not rows:
        raise ValueError(
            "No variantEvaluations were found."
        )
    df = pl.from_dicts(
        rows,
        infer_schema_length=None
    )
    # ========================================================
    # Canonical columns first
    # ========================================================
    columns = [
        "CONTIG",
        "START",
        "REF",
        "ALT",
        "EXOMISER_JSON_GENES",
    ]
    # ========================================================
    # Frequency columns
    #
    # AF only
    # ========================================================
    for source in frequency_sources:
        column = (
            f"{source}_AF"
        )
        if column in df.columns:
            columns.append(
                column
            )
    # ========================================================
    # Pathogenicity predictor columns
    # ========================================================
    for source in pathogenicity_sources:
        raw_column = (
            f"{source}_RAW_SCORE"
        )
        scaled_column = (
            f"{source}_SCALED_SCORE"
        )
        score_column = (
            f"{source}_SCORE"
        )
        if raw_column in df.columns:
            columns.append(
                raw_column
            )
        if scaled_column in df.columns:
            columns.append(
                scaled_column
            )
        if score_column in df.columns:
            columns.append(
                score_column
            )
    # ========================================================
    # Exomiser supporting fields
    # ========================================================
    supporting_columns = [
        "FREQUENCY_SCORE",
        "PATHOGENICITY_SCORE",
        "PREDICTED_PATHOGENIC",
        "PASSED_FILTER_TYPES",
        "COMPATIBLE_INHERITANCE_MODES",
        "CONTRIBUTING_INHERITANCE_MODES",
    ]
    for column in supporting_columns:
        if column in df.columns:
            columns.append(
                column
            )
    # ========================================================
    # Keep any unexpected useful fields at the end
    # ========================================================
    remaining_columns = [
        column
        for column in df.columns
        if column not in columns
    ]
    df = df.select(
        columns
        +
        remaining_columns
    )
    return df

# ============================================================
# Main
# ============================================================
def main():
    args = parse_arguments()
    input_path = Path(args.input)
    output_path = Path(args.output)

    try:
        # ----------------------------------------------------
        # Validate input
        # ----------------------------------------------------
        if not input_path.is_file():
            raise FileNotFoundError(
                f"Input file not found: "
                f"{input_path}"
            )
        output_path.parent.mkdir(
            parents=True,
            exist_ok=True
        )
        
        # ----------------------------------------------------
        # Read JSON
        # ----------------------------------------------------
        gene_records = (
            read_exomiser_json(
                input_path
            )
        )
        
        # ----------------------------------------------------
        # Extract variants
        # ----------------------------------------------------
        (
            rows,
            frequency_sources,
            pathogenicity_sources,
        ) = extract_variants(
            gene_records
        )

        # ----------------------------------------------------
        # Build dataframe
        # ----------------------------------------------------
        df = create_dataframe(
            rows,
            frequency_sources,
            pathogenicity_sources,
        )

        # ----------------------------------------------------
        # Write TSV
        # ----------------------------------------------------
        df.write_csv(
            output_path,
            separator="\t"
        )

        # ----------------------------------------------------
        # Summary
        # ----------------------------------------------------
        print(
            f"Gene records          : "
            f"{len(gene_records)}"
        )
        print(
            f"Unique variants       : "
            f"{df.height}"
        )
        print(
            "Frequency sources      : "
            + (
                ", ".join(
                    frequency_sources
                )
                if frequency_sources
                else "None"
            )
        )
        print(
            "Pathogenicity sources  : "
            + (
                ", ".join(
                    pathogenicity_sources
                )
                if pathogenicity_sources
                else "None"
            )
        )
        print(
            f"Output                : "
            f"{output_path}"
        )
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