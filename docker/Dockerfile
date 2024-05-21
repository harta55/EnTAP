FROM --platform=linux/amd64 python:3.10
ENV DEBIAN_FRONTEND=noninteractive

# Update and install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        wget \
        unzip \
        gnupg \
        default-jre \
        python3-pip \
        cmake \
        build-essential \
        zlib1g-dev \
        r-base \
        elfutils \
        libtclap-dev \
        libgsl-dev \
        sqlite3 \
        python3-numpy \
        python3-pandas \
        python3-xmltodict \
        python3-matplotlib \
        python3-setuptools && \
        apt-get clean && \
        rm -rf /var/lib/apt/lists/*

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

# Download TransDecoder dependencies
RUN curl -L https://cpanmin.us | perl - App::cpanminus \
&& cpanm install DB_File \
&& cpanm install URI::Escape

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
&& make install \
&& cd .. \
&& rm -rf RSEM-1.3.3 v1.3.3.tar.gz

# Download eggnog-mapper
RUN wget https://github.com/eggnogdb/eggnog-mapper/archive/refs/tags/2.1.12.tar.gz \
    && tar -zxvf 2.1.12.tar.gz \
    && cd eggnog-mapper-2.1.12 \
    && python setup.py install \
    && cd .. \
    && rm -rf eggnog-mapper-2.1.12 2.1.12.tar.gz

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
