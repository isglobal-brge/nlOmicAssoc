#################################################################
# Dockerfile
#
# Version:          1
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
    libcurl4-openssl-dev 
    
    
#Rcompression needed by RCurl
RUN install2.r -r "http://www.omegahat.net/R" --deps TRUE \
    Rcompression \
    && rm -rf /tmp/downloaded_packages/

## Finally ready to install the R packages.  NOTE: failure to install a package doesn't throw an image build error.
RUN install2.r --error --deps TRUE \
    RCurl \
    mfp \
    stringr \
    RUnit \
    && rm -rf /tmp/downloaded_packages/

## Add biocLite to install Biobase
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("Biobase");'


##That's all for the moment
