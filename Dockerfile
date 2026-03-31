FROM rocker/r-ver:4.4.0

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV R_LIBS_USER=/app/.r_libs

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    gfortran \
    git \
    curl \
    ca-certificates \
    locales \
    pkg-config \
    python3 \
    python3-pip \
    libblas-dev \
    liblapack-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfftw3-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgdal-dev \
    libgeos-dev \
    libgit2-dev \
    libglpk-dev \
    libgmp3-dev \
    libgsl-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    libjpeg-dev \
    liblzma-dev \
    libpcre2-dev \
    libpng-dev \
    libproj-dev \
    libsqlite3-dev \
    libssl-dev \
    libtiff5-dev \
    libudunits2-dev \
    libxml2-dev \
    libxt-dev \
    zlib1g-dev \
 && sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen \
 && locale-gen \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt /app/requirements.txt
COPY dependencies /app/dependencies
COPY R/install_exact_dependencies.R /app/R/install_exact_dependencies.R

RUN python3 -m pip install --no-cache-dir -r /app/requirements.txt
RUN mkdir -p /app/.r_libs
RUN Rscript --vanilla /app/R/install_exact_dependencies.R

COPY . /app

ENTRYPOINT ["python3", "main_pipeline.py"]
CMD ["--help"]
