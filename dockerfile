# ...
FROM mambaorg/micromamba:1.5.8
ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]

ENV DEBIAN_FRONTEND=noninteractive LC_ALL=C.UTF-8 LANG=C.UTF-8 TZ=UTC

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl ca-certificates unzip procps pigz && \
    rm -rf /var/lib/apt/lists/*

# 1) Ensure strict channel priority (reduces conflicts)
RUN micromamba config set channel_priority strict && \
    micromamba config append channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config append channels defaults

# 2) Create env with compatible versions (Java included for fastqc/bbmap)
RUN micromamba create -y -n bio \
      python=3.10 \
      openjdk=17 \
      fastp=0.23.4 \
      fastqc=0.12.1 \
      multiqc=1.23 \
      seqkit=2.8.2 \
      bbmap=39.01 && \
    micromamba clean -a -y

# 3) Install AWS CLI v2 (outside conda)
RUN curl -fsSL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip" && \
    unzip -q /tmp/awscliv2.zip -d /tmp && \
    /tmp/aws/install && rm -rf /tmp/aws /tmp/awscliv2.zip

# 4) PATH + (optional) make /usr/bin/java to silence JAVA_CMD warnings
ENV PATH="/opt/conda/envs/bio/bin:/usr/local/bin:${PATH}"
RUN ln -sf /opt/conda/envs/bio/bin/java /usr/bin/java || true

# Non-root user + workspace
ARG USERNAME=nfuser UID=1000
RUN useradd -m -u ${UID} -s /bin/bash ${USERNAME}
USER ${USERNAME}
WORKDIR /workspace