# Base image.
FROM ubuntu:18.04

# Set noninteractive mode for apt-get.
ENV DEBIAN_FRONTEND noninteractive


# Update and install basic packages.
RUN apt-get update -qq \
  && apt-get install -qq -y curl git unzip wget gnupg


# Install Oracle Java 8 which is used by InterProScan. Rather than use
# the apt-add-repository program (which installs a lot of packages) we will
# manually add the repository.
RUN echo "deb http://ppa.launchpad.net/openjdk-r/ppa/ubuntu xenial main" >> /etc/apt/sources.list \
  && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 86F44E2A

RUN apt-get update \
  && apt-get -y install openjdk-8-jdk \
  && update-alternatives --config java

# Install Python3
RUN apt-get install -qq -y python3-pip \
&& pip3 install -q numpy pandas xmltodict matplotlib \
&& ln -s /usr/bin/python3 /usr/bin/python

# Install R
RUN apt-get -y install software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
&& add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/' \
&& apt-get update \
&& apt-get -y install r-base \
&& echo "install.packages('ggplot2',repos='https://ftp.osuosl.org/pub/cran/')" | R --vanilla \
&& echo "install.packages('BiocManager',repos='https://ftp.osuosl.org/pub/cran/')" | R --vanilla \
&& echo "BiocManager::install('seqLogo')" | R --vanilla

# TransDecoder
RUN git clone https://github.com/TransDecoder/TransDecoder.git \
&& cd TransDecoder \
&& git checkout TransDecoder-v5.3.0 \
&& cp -R . /usr/lib/TransDecoder-v5.3.0

# Install InterProScan
RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.36-75.0/alt/interproscan-core-5.36-75.0.tar.gz \
  && tar -zxvf interproscan-core-5.36-75.0.tar.gz \
  && mv interproscan-5.36-75.0 /usr/local/interproscan

# Install IPRScan dependency
RUN apt-get install -y elfutils

# Install Diamond
RUN wget https://github.com/bbuchfink/diamond/releases/download/v0.9.19/diamond-linux64.tar.gz \
  && tar -zxvf diamond-linux64.tar.gz \
  && mv diamond /usr/local/bin/diamond

# Install RSEM
RUN apt-get -y install zlib1g-dev
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz \
  && tar -zxvf v1.3.0.tar.gz \
  && cd RSEM-1.3.0 \
  && make \
  && make ebseq \
  && make install

# Install EnTAP
RUN apt-get -y install cmake
COPY . /tmp/entap
RUN cd /tmp/entap \
  && cmake CMakeLists.txt \
  && make \
  && make install \
  && cp src/entap_graphing.py /usr/local/bin

COPY entap_config.ini /entap_config.ini


# Tini for signal processing and zombie killing.
ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini


# Cleanup
RUN rm -rf RSEM-1.3.0
RUN rm -fr TransDecoder
RUN rm diamond-linux64.tar.gz && rm diamond_manual.pdf
RUN rm interproscan-core-5.36-75.0.tar.gz
RUN rm v1.3.0.tar.gz
RUN rm LICENSE
RUN rm -fr /tmp/entap

# Define the command and parameters that will be executed when this
# container is first run.
ENTRYPOINT ["/tini", "--"]
