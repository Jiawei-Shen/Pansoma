# syntax=docker/dockerfile:1

FROM nvidia/cuda:12.1.1-cudnn8-devel-ubuntu22.04

ARG VG_VERSION=1.74.1
ARG TARGETARCH

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONNOUSERSITE=1 \
    PANSOMA_HOME=/opt/pansoma \
    PANSOMA_PYTHON=/opt/venv/bin/python \
    VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:$PATH

RUN if [ -n "${TARGETARCH}" ] && [ "${TARGETARCH}" != "amd64" ]; then \
        echo "Pansoma's bundled vg release binary requires linux/amd64" >&2; \
        exit 1; \
    fi

RUN apt-get update && apt-get install -y --no-install-recommends \
        bcftools \
        build-essential \
        ca-certificates \
        curl \
        jq \
        libgomp1 \
        python3 \
        python3-dev \
        python3-pip \
        python3-venv \
        samtools \
        tabix \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv "${VIRTUAL_ENV}" \
    && "${VIRTUAL_ENV}/bin/python" -m pip install --no-cache-dir \
        --upgrade pip setuptools wheel

RUN curl --fail --location --retry 3 \
        "https://github.com/vgteam/vg/releases/download/v${VG_VERSION}/vg" \
        --output /usr/local/bin/vg \
    && chmod 0755 /usr/local/bin/vg \
    && vg version

COPY docker/requirements.pipeline.txt /tmp/requirements.pipeline.txt
RUN "${PANSOMA_PYTHON}" -m pip install --no-cache-dir \
        --extra-index-url https://download.pytorch.org/whl/cu121 \
        torch==2.1.2+cu121 torchvision==0.16.2+cu121 \
    && "${PANSOMA_PYTHON}" -m pip install --no-cache-dir \
        -r /tmp/requirements.pipeline.txt

WORKDIR /opt/pansoma
COPY . /opt/pansoma

RUN "${PANSOMA_PYTHON}" -m pip install \
        --no-cache-dir --no-deps --no-build-isolation --editable . \
    && PANSOMA_PYTHON="${PANSOMA_PYTHON}" bash scripts/build_fast_writer.sh \
    && "${PANSOMA_PYTHON}" scripts/check_environment.py --strict-external --ml

RUN install -m 0755 docker/entrypoint.sh /usr/local/bin/pansoma-entrypoint \
    && ln -s /usr/local/bin/pansoma-entrypoint /usr/local/bin/pansoma \
    && mkdir -p /work /output

ENTRYPOINT ["/usr/local/bin/pansoma-entrypoint"]
CMD ["--help"]
