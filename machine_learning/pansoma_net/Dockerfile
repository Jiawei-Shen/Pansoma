# Use NVIDIA CUDA base image
FROM nvidia/cuda:12.1.0-cudnn8-devel-ubuntu20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    PATH="/opt/conda/bin:$PATH"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    curl \
    ca-certificates \
    libgl1-mesa-glx \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
WORKDIR /opt
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Set the default shell to bash for Conda
SHELL ["/bin/bash", "-c"]

# Create and activate a minimal Conda environment with Python 3.8
RUN conda create -n pansoma_net_env python=3.8 -y && \
    conda clean --all -y
    
# Set workspace as working directory
WORKDIR /workspace

# Copy project files into the container
# COPY . .
COPY requirements.txt /workspace/requirements.txt

# Activate Conda environment and install dependencies via pip
#COPY requirements.txt /workspace/requirements.txt
RUN conda run -n pansoma_net_env pip install --no-cache-dir -r /workspace/requirements.txt

# Set the container entrypoint to activate Conda environment
ENTRYPOINT ["/bin/bash", "-c", "cd /workspace/ && source /opt/conda/bin/activate pansoma_net_env && bash"]
# CMD ["/bin/bash"]


