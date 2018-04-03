#################################################################
# Dockerfile
#
# Version:          1
# Date              160617 
# Review1           030418 added caret and gam removed mda from r-base: 3.4.3 
# Software:         R
# Description:      R and necessary packages 
# Website:          https://github.com/isglobal-brge/nlOmicAssoc|https://hub.docker.com/r/lnonell/nlomicsassoc
# Tags:             None, for the moment
# Base Image:       R
#################################################################

##Image created on a debian
FROM r-base:3.4.0

#xml needed by Rcompression
#curl needed by RCurl
RUN apt-get update && apt-get install -y \
    r-cran-xml \
    #curl \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev #per paquet xml2

ENV PATH=pkg-config:$PATH   

#Rcompression needed by RCurl
RUN install2.r -r http://www.omegahat.net/R --deps TRUE \
    Rcompression \
    && rm -rf /tmp/downloaded_packages/
    
#splines has been moved to archive
RUN wget https://cran.r-project.org/src/contrib/Archive/splines/splines_2.0-7.tar.gz -P /tmp/downloaded_packages
    

## Finally ready to install the R packages.  NOTE: failure to install a package doesn't throw an image build error.
RUN install2.r --error --deps TRUE \
    RCurl \
    xml2 \
    mfp \
    stringr \
    RUnit \
    mgcv \
    mvtnorm \

    && rm -rf /tmp/downloaded_packages/
    
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Biobase");biocLite("KEGGgraph");biocLite("biomaRt");'
    
## Finally ready to install the R packages.  NOTE: failure to install a package doesn't throw an image build error.
RUN install2.r --error --deps TRUE \
    MASS \
    gam \
    mboost \
    partDSA \
    nnet \
    randomForest \
    vegan \
    devtools \
    NeuralNetTools \
    mice \
    caret \
    && rm -rf /tmp/downloaded_packages/


## Add biocLite to install Biobase
RUN Rscript -e 'install.packages("/tmp/downloaded_packages/splines_2.0-7.tar.gz", repos = NULL, type = "source")'


##That's all for the moment

