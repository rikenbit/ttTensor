# Base Image
FROM bioconductor/bioconductor_docker:devel

# Install R Packages
RUN R -e "devtools::install_github('rikenbit/ttTensor', \
    upgrade='always', force=TRUE, INSTALL_opts = '--install-tests');\
    tools::testInstalledPackage('ttTensor')"