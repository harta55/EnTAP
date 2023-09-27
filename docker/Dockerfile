FROM --platform=linux/amd64 debian:latest

ENV DEBIAN_FRONTEND=noninteractive

# Update and install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        wget \
        unzip \
        gnupg \
		git \
        default-jre \
        python3-pip \
        cmake \
        build-essential \
        zlib1g-dev \
        r-base \
        elfutils \
        libtclap-dev \
        libgsl-dev \
        python3-numpy \
        python3-pandas \
        python3-xmltodict \
        python3-matplotlib

# Create a symbolic link for Python
RUN ln -s /usr/bin/python3 /usr/bin/python && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install InterProScan
RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.64-96.0/alt/interproscan-core-5.64-96.0.tar.gz \
&& tar -zxvf interproscan-core-5.64-96.0.tar.gz \
&& mv interproscan-5.64-96.0 /usr/local/interproscan

# Install IPRScan dependency
RUN apt-get install -y elfutils

# Download Diamond
RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz && \
    tar -zxvf diamond-linux64.tar.gz && \
    mv diamond /usr/local/bin/diamond && \
    rm diamond-linux64.tar.gz

# Download TransDecoder
RUN git clone https://github.com/TransDecoder/TransDecoder.git \
&& cd TransDecoder \
&& git checkout TransDecoder-v5.7.1 \
&& cp -R . /usr/lib/TransDecoder-v5.7.1

# Install R and required packages
RUN R -e "install.packages('ggplot2', repos='https://ftp.osuosl.org/pub/cran/')" && \
    R -e "install.packages('BiocManager', repos='https://ftp.osuosl.org/pub/cran/')" && \
    R -e "BiocManager::install('seqLogo')"

# Install RSEM
RUN apt-get -y install zlib1g-dev \
&& wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz \
&& tar -zxvf v1.3.3.tar.gz \
&& cd RSEM-1.3.3 \
&& make \
&& make ebseq \
&& make install

# Clone EnTAP repository and install requirements
RUN git clone https://github.com/harta55/EnTAP.git /tmp/entap && \
    cd /tmp/entap && \
    cmake CMakeLists.txt && \
    make && \
    make install && \
	mv src/entap_graphing.py /usr/local/bin/entap_graphing.py && \
    cd .. && \
    rm -rf /tmp/entap

CMD ["bash"]