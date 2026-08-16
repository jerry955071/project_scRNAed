FROM rocker/rstudio:4.4.3

USER root

ARG DEBIAN_FRONTEND=noninteractive

# System libraries required by Seurat and the statistical packages.
# cmake + libnlopt-dev are especially important for nloptr/lme4/car.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        cmake \
        pkg-config \
        pandoc \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libnlopt-dev \
        libcairo2-dev \
        libxt-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        libjpeg-dev \
        libtiff-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libglpk-dev \
        libgsl-dev \
        libicu-dev \
    && rm -rf /var/lib/apt/lists/*


# Install Matrix first.
#
# This is intentional:
# TMB, glmmTMB, lme4, SeuratObject, etc. should all be installed
# against the same Matrix version from Rocker's CRAN snapshot.
RUN install2.r \
        --error \
        -n 4 \
        Matrix


# Packages required by run-REDIT-regression.Rmd
# plus Seurat.
RUN install2.r \
        --error \
        --skipinstalled \
        -n 4 \
        Seurat \
        glmmTMB \
        doParallel \
        foreach \
        DHARMa \
        emmeans \
        car \
        dplyr \
        magrittr \
        ggplot2 \
        patchwork \
        ComplexUpset \
        data.table \
        scales \
        knitr \
        rmarkdown


# Build-time validation.
RUN Rscript -e ' \
    pkgs <- c( \
        "Seurat", \
        "Matrix", \
        "glmmTMB", \
        "TMB", \
        "lme4", \
        "doParallel", \
        "foreach", \
        "DHARMa", \
        "emmeans", \
        "car", \
        "dplyr", \
        "magrittr", \
        "ggplot2", \
        "patchwork", \
        "ComplexUpset", \
        "data.table", \
        "scales", \
        "knitr", \
        "rmarkdown" \
    ); \
    failed <- pkgs[!vapply( \
        pkgs, \
        requireNamespace, \
        FUN.VALUE = logical(1), \
        quietly = TRUE \
    )]; \
    if (length(failed)) { \
        stop("Missing packages: ", paste(failed, collapse = ", ")); \
    }; \
    if (!rmarkdown::pandoc_available()) { \
        stop("Pandoc is not available"); \
    }; \
    cat("\nInstalled package versions:\n"); \
    print(vapply( \
        pkgs, \
        function(x) as.character(packageVersion(x)), \
        FUN.VALUE = character(1) \
    )); \
    cat("\nR version:\n"); \
    print(R.version.string); \
    cat("\nCRAN repository:\n"); \
    print(getOption("repos")) \
'


EXPOSE 8787/tcp

CMD ["/init"]