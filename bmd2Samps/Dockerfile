FROM bioconductor/bioconductor_docker
# To setup other dependencies

#COPY setup.sh /setup.sh
RUN Rscript -e "install.packages('argparse')"

COPY . bmd2Samps
WORKDIR bmd2Samps

ENTRYPOINT ["Rscript","buildv1database.R"]
VOLUME ["/tmp"]
