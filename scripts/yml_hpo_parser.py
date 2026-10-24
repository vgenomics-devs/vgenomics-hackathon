#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

HPO_PATTERN = re.compile(r"^HP:\d{7}$")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Update hpoIds in an Exomiser YAML configuration file."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input Exomiser YAML template file."
    )

    parser.add_argument(
        "-p",
        "--hpo-ids",
        required=True,
        help=(
            "Comma-separated HPO IDs. "
            "Example: HP:0001250,HP:0001263,HP:0004322"
        )
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output YAML file."
    )

    return parser.parse_args()

def validate_hpo_ids(hpo_string):
    raw_ids = [
        value.strip()
        for value in hpo_string.split(",")
        if value.strip()
    ]

    if not raw_ids:
        raise ValueError("No HPO IDs were provided.")

    valid_ids = []
    invalid_ids = []

    for hpo_id in raw_ids:
        if HPO_PATTERN.fullmatch(hpo_id):
            valid_ids.append(hpo_id)
        else:
            invalid_ids.append(hpo_id)

    if invalid_ids:
        raise ValueError(
            "Invalid HPO ID(s): "
            + ", ".join(invalid_ids)
            + ". Expected format: HP:0000000"
        )

    # Remove duplicates while preserving original order
    valid_ids = list(dict.fromkeys(valid_ids))

    return valid_ids


def update_hpo_ids(input_file, output_file, hpo_ids):
    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.exists():
        raise FileNotFoundError(
            f"Input YAML file not found: {input_path}"
        )

    if not input_path.is_file():
        raise ValueError(
            f"Input path is not a file: {input_path}"
        )

    with input_path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    formatted_hpo_ids = ", ".join(hpo_ids)

    updated = False

    for index, line in enumerate(lines):
        stripped = line.lstrip()

        if stripped.startswith("hpoIds:"):
            indentation = line[: len(line) - len(stripped)]

            lines[index] = (
                f"{indentation}hpoIds: [{formatted_hpo_ids}]\n"
            )

            updated = True
            break

    if not updated:
        raise ValueError(
            "Could not find an 'hpoIds:' entry in the input YAML file."
        )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with output_path.open("w", encoding="utf-8") as handle:
        handle.writelines(lines)


def main():
    args = parse_arguments()

    try:
        hpo_ids = validate_hpo_ids(
            args.hpo_ids
        )

        update_hpo_ids(
            input_file=args.input,
            output_file=args.output,
            hpo_ids=hpo_ids
        )

        print(
            f"Updated {len(hpo_ids)} HPO IDs successfully."
        )

        print(
            f"Output file: {args.output}"
        )

        return 0

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr
        )

        return 1


if __name__ == "__main__":
    sys.exit(main())
