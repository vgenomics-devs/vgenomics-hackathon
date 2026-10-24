#!/usr/bin/env python3
import argparse
import gzip
import sys
from pathlib import Path

import polars as pl

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Extract selected fields from a normalized/split single-sample "
            "VCF.gz file and write them to TSV."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input normalized/split VCF.gz file."
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file."
    )

    return parser.parse_args()


def open_vcf(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "r", encoding="utf-8")


def get_format_value(format_string, sample_string, key):

    if not format_string or not sample_string:
        return None

    if format_string == "." or sample_string == ".":
        return None

    keys = format_string.split(":")
    values = sample_string.split(":")

    try:
        idx = keys.index(key)
    except ValueError:
        return None

    if idx >= len(values):
        return None

    value = values[idx]

    if value in ("", "."):
        return None

    return value


def get_zygosity(gt):

    if gt is None or gt in ("", ".", "./.", ".|."):
        return "UNKNOWN"

    # Treat phased and unphased genotypes equivalently
    # for zygosity classification.
    alleles = gt.replace("|", "/").split("/")

    if any(allele == "." for allele in alleles):
        return "UNKNOWN"

    # Haploid genotype
    # 0 -> reference allele
    # 1,2,... -> alternate allele
    if len(alleles) == 1:

        if alleles[0] == "0":
            return "HAPLOID_REF"

        return "HAPLOID_ALT"

    unique_alleles = set(alleles)

    # 0/0
    if unique_alleles == {"0"}:
        return "HOM_REF"

    # 1/1, 2/2, etc.
    if len(unique_alleles) == 1 and "0" not in unique_alleles:
        return "HOM_ALT"

    # 0/1, 0/2, etc.
    if "0" in unique_alleles and len(unique_alleles) > 1:
        return "HET"

    # 1/2, 1/3, etc.
    if "0" not in unique_alleles and len(unique_alleles) > 1:
        return "HET_ALT"

    return "UNKNOWN"


def calculate_vaf(ad):

    if ad is None or ad in ("", "."):
        return None

    try:
        depths = [
            int(value) if value not in ("", ".") else 0
            for value in ad.split(",")
        ]

    except ValueError:
        return None

    if len(depths) < 2:
        return None

    total_depth = sum(depths)

    if total_depth == 0:
        return None

    # Expected after:
    # bcftools norm -m -any
    #
    # AD = REF_DEPTH,ALT_DEPTH
    alt_depth = depths[1]

    return alt_depth / total_depth


def extract_vcf(input_vcf):

    records = []
    sample_name = None
    header_found = False

    with open_vcf(input_vcf) as vcf:

        for line in vcf:

            # Skip VCF metadata
            if line.startswith("##"):
                continue

            # VCF column header
            if line.startswith("#CHROM"):

                header = line.rstrip("\n").split("\t")

                if len(header) < 10:
                    raise ValueError(
                        "VCF does not contain a sample column."
                    )

                sample_name = header[9]
                header_found = True

                continue

            if line.startswith("#"):
                continue

            if not header_found:
                raise ValueError(
                    "VCF #CHROM header was not found before variant records."
                )

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                raise ValueError(
                    f"Malformed VCF record with fewer than 10 columns: {line.strip()}"
                )

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_value = fields[6]
            info = fields[7]

            format_string = fields[8]
            sample_string = fields[9]

            gt = get_format_value(
                format_string,
                sample_string,
                "GT"
            )

            ad = get_format_value(
                format_string,
                sample_string,
                "AD"
            )

            dp = get_format_value(
                format_string,
                sample_string,
                "DP"
            )

            gq = get_format_value(
                format_string,
                sample_string,
                "GQ"
            )

            pl = get_format_value(
                format_string,
                sample_string,
                "PL"
            )

            zygosity = get_zygosity(gt)
            vaf = calculate_vaf(ad)

            records.append(
                {
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "QUAL": qual,
                    "FILTER": filter_value,
                    "INFO": info,

                    "GT": gt,
                    "ZYGOSITY": zygosity,

                    "AD": ad,
                    "VAF": vaf,

                    "DP": dp,
                    "GQ": gq,
                    "PL": pl,
                }
            )

    if not header_found:
        raise ValueError(
            "VCF header line '#CHROM' was not found."
        )

    return records, sample_name


def create_dataframe(records):

    if not records:
        return pl.DataFrame(
            schema={
                "CHROM": pl.String,
                "POS": pl.Int64,
                "REF": pl.String,
                "ALT": pl.String,
                "QUAL": pl.Float64,
                "FILTER": pl.String,
                "INFO": pl.String,
                "GT": pl.String,
                "ZYGOSITY": pl.String,
                "AD": pl.String,
                "VAF": pl.Float64,
                "DP": pl.Int64,
                "GQ": pl.Float64,
                "PL": pl.String,
            }
        )

    df = pl.DataFrame(records)

    df = df.with_columns(

        pl.col("POS")
        .cast(pl.Int64, strict=False),

        pl.col("QUAL")
        .replace(".", None)
        .cast(pl.Float64, strict=False),

        pl.col("DP")
        .replace(".", None)
        .cast(pl.Int64, strict=False),

        pl.col("GQ")
        .replace(".", None)
        .cast(pl.Float64, strict=False),

        pl.col("VAF")
        .cast(pl.Float64, strict=False),

    )

    return df


def main():

    args = parse_arguments()

    input_path = Path(args.input)
    output_path = Path(args.output)

    try:

        # ----------------------------------------------------
        # Input validation
        # ----------------------------------------------------

        if not input_path.exists():
            raise FileNotFoundError(
                f"Input VCF not found: {input_path}"
            )

        if not input_path.is_file():
            raise ValueError(
                f"Input path is not a file: {input_path}"
            )

        if not str(input_path).endswith(".vcf.gz"):
            raise ValueError(
                "Input file must be a gzipped VCF ending with '.vcf.gz'."
            )

        # ----------------------------------------------------
        # Create output directory if required
        # ----------------------------------------------------

        output_path.parent.mkdir(
            parents=True,
            exist_ok=True
        )

        # ----------------------------------------------------
        # Extract VCF
        # ----------------------------------------------------

        records, sample_name = extract_vcf(
            input_path
        )

        # ----------------------------------------------------
        # Convert to Polars DataFrame
        # ----------------------------------------------------

        df = create_dataframe(
            records
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

        print(f"Sample   : {sample_name}")
        print(f"Variants : {df.height}")
        print(f"Output   : {output_path}")

        return 0

    except Exception as error:

        print(
            f"ERROR: {error}",
            file=sys.stderr
        )

        return 1

if __name__ == "__main__":
    sys.exit(main())
