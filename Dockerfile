FROM rocker/r2u:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV VENV_PATH=/opt/customgwas-venv

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    python3-dev \
    build-essential \
    git \
    curl \
    ca-certificates \
    openjdk-17-jre-headless \
    libopenblas-dev \
    liblapack-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    bcftools \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK with distro-dependent package names.
RUN apt-get update && \
    (apt-get install -y --no-install-recommends plink || apt-get install -y --no-install-recommends plink2) && \
    rm -rf /var/lib/apt/lists/*

# Normalize PLINK command name for notebook checks (which('plink')).
RUN if ! command -v plink >/dev/null 2>&1 && command -v plink2 >/dev/null 2>&1; then \
      ln -sf "$(command -v plink2)" /usr/local/bin/plink; \
    fi

COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv "$VENV_PATH" \
    && "$VENV_PATH/bin/python" -m pip install --no-cache-dir --upgrade pip setuptools wheel \
    && "$VENV_PATH/bin/python" -m pip install --no-cache-dir -r /tmp/requirements.txt

ENV PATH="$VENV_PATH/bin:$PATH"

COPY docker/install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R

COPY beagle.27Feb25.75f.jar /work/beagle.27Feb25.75f.jar

COPY docker/entrypoint.sh /usr/local/bin/customgwas-entrypoint.sh
RUN chmod +x /usr/local/bin/customgwas-entrypoint.sh

WORKDIR /work

EXPOSE 8888
ENTRYPOINT ["/usr/local/bin/customgwas-entrypoint.sh"]

