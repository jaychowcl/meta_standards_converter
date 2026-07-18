# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
ARG DOCKER_CLI_VERSION=29.6.2
FROM docker:${DOCKER_CLI_VERSION}-cli AS docker_cli

FROM python:3.12-slim

ARG NEXTFLOW_VERSION=26.04.2
ARG NEXTFLOW_SHA256=52a2ce22be15d747369a70050339443fc325005bea59a49e0ee25d96fae9cc51

COPY --from=docker_cli /usr/local/bin/docker /usr/local/bin/docker

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        bash \
        ca-certificates \
        curl \
        git \
        gffread \
        openjdk-21-jre-headless \
    && curl --fail --location --show-error \
        --output /tmp/nextflow \
        "https://github.com/nextflow-io/nextflow/releases/download/v${NEXTFLOW_VERSION}/nextflow-${NEXTFLOW_VERSION}-dist" \
    && echo "${NEXTFLOW_SHA256}  /tmp/nextflow" | sha256sum --check --strict \
    && install --mode=0755 /tmp/nextflow /usr/local/bin/nextflow \
    && rm --force /tmp/nextflow \
    && rm --recursive --force /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml README.md LICENSE ./
COPY src ./src

RUN pip install --no-cache-dir '.[h5ad]'

RUN java -version \
    && nextflow -version \
    && docker --version

CMD ["geo2ae", "--help"]
