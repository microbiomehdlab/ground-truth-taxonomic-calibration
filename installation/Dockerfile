# --- Base (pin micromamba image) ---
FROM mambaorg/micromamba:1.5.6
SHELL ["/bin/bash", "-lc"]
WORKDIR /opt

# --- Create env with pinned tool versions ---
# Channel order: conda-forge first (runtime libs), then bioconda (bio tools)
RUN micromamba create -y -n taxonomic_tools -c conda-forge -c bioconda \
      python=3.11 pip \
      kraken2=2.1.6 \
      bracken=3.1 \
      metaphlan=4.2.2 \
      bowtie2=2.5.4 \
      art=2016.06.05 \
      seqtk=1.4 \
      fastqc=0.12.1 \
      fastp=0.23.4 \
  && micromamba clean --all --yes

# --- Runtime certificates ---
USER root
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# --- Back to non-root micromamba user ---
ARG MAMBA_USER=mambauser
USER $MAMBA_USER

# Optional: convenient activation when you open a shell
ENV MAMBA_DOCKERFILE_ACTIVATE=1
RUN echo 'eval "$(micromamba shell hook --shell bash)"' >> ~/.bashrc && \
    echo 'micromamba activate taxonomic_tools' >> ~/.bashrc

# --- Always execute inside the conda env ---
# This lets you run tools directly, e.g.:
#   docker run IMAGE kraken2 --version
ENTRYPOINT ["micromamba","run","-n","taxonomic_tools"]
CMD ["bash"]
