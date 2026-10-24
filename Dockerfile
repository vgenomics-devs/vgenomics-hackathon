FROM ubuntu:24.04

# ============================================================
# BUILD ARGUMENTS
# ============================================================
ARG DEBIAN_FRONTEND=noninteractive
ARG EXOMISER_VERSION=15.1.0

# ============================================================
# SYSTEM PACKAGES
# ============================================================
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        openjdk-21-jre-headless \
        python3 \
        python3-pip \
        bcftools \
        tabix \
        wget \
        curl \
        unzip \
        ca-certificates \
        bash \
        coreutils \
        procps \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# PYTHON PACKAGES
# ============================================================
RUN pip3 install \
    --no-cache-dir \
    --break-system-packages \
    polars

# ============================================================
# INSTALL EXOMISER 15.1.0
# ============================================================
WORKDIR /opt
RUN wget \
    -O /tmp/exomiser.zip \
    "https://github.com/exomiser/Exomiser/releases/download/${EXOMISER_VERSION}/exomiser-cli-${EXOMISER_VERSION}-distribution.zip" \
    && unzip /tmp/exomiser.zip -d /opt \
    && rm /tmp/exomiser.zip \
    && mv "/opt/exomiser-cli-${EXOMISER_VERSION}" /opt/exomiser

# ============================================================
# CREATE PIPELINE DIRECTORIES
# ============================================================
RUN mkdir -p \
    /pipeline/scripts \
    /pipeline/config \
    /input \
    /output \
    /reference \
    /exomiser-data

# ============================================================
# COPY PIPELINE SCRIPTS
# ============================================================
COPY scripts/yml_hpo_parser.py \
    /pipeline/scripts/yml_hpo_parser.py

COPY scripts/wgs_vcf_extract.py \
    /pipeline/scripts/wgs_vcf_extract.py

COPY scripts/exomiser_json_extract.py \
    /pipeline/scripts/exomiser_json_extract.py

COPY scripts/exomiser_json_tsv_merge.py \
    /pipeline/scripts/exomiser_json_tsv_merge.py

COPY scripts/wgs_vcf_exomiser_merge.py \
    /pipeline/scripts/wgs_vcf_exomiser_merge.py

COPY scripts/wgs_variants_prioritizer.py \
    /pipeline/scripts/wgs_variants_prioritizer2.py

# ============================================================
# COPY MAIN PIPELINE SCRIPT
# ============================================================
COPY run_wgs_pipeline.sh \
    /pipeline/run_wgs_pipeline.sh

# ============================================================
# PERMISSIONS
# ============================================================
RUN chmod +x \
    /pipeline/run_wgs_pipeline.sh \
    /pipeline/scripts/*.py

# ============================================================
# VERIFY INSTALLATIONS DURING BUILD
# ============================================================
RUN java -version && \
    python3 --version && \
    bcftools --version | head -1 && \
    tabix --version | head -1 && \
    python3 -c "import polars as pl; print('Polars:', pl.__version__)" && \
    test -f /opt/exomiser/exomiser-cli-${EXOMISER_VERSION}.jar

# ============================================================
# WORKING DIRECTORY
# ============================================================
WORKDIR /pipeline

# ============================================================
# ENTRYPOINT
# ============================================================
ENTRYPOINT ["/pipeline/run_wgs_pipeline.sh"]