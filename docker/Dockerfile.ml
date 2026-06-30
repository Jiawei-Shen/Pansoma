# Pansoma ML container.
#
# Build from the repository root:
#   docker build -f docker/Dockerfile.ml -t pansoma-ml .

FROM nvidia/cuda:12.1.0-cudnn8-devel-ubuntu20.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    PATH="/opt/conda/bin:$PATH"

RUN apt-get update && apt-get install -y \
    wget \
    git \
    curl \
    ca-certificates \
    libgl1-mesa-glx \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

SHELL ["/bin/bash", "-c"]

RUN conda create -n pansoma_ml python=3.8 -y && \
    conda clean --all -y

WORKDIR /workspace/Pansoma

COPY ["machine_learning/pansoma_net/requirements.txt", "/tmp/pansoma_ml_requirements.txt"]
RUN conda run -n pansoma_ml pip install --no-cache-dir -r /tmp/pansoma_ml_requirements.txt

ENTRYPOINT ["/bin/bash", "-c", "source /opt/conda/bin/activate pansoma_ml && cd /workspace/Pansoma && bash"]

