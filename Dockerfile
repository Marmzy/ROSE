FROM python:3.10.13-slim

WORKDIR /ROSE

COPY ROSE.py requirements.txt README.md LICENSE ./

RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
    bzip2 \
    gcc \
    git \
    make \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \ 
    libncursesw5-dev \
    wget \
    zlib1g-dev \
    && apt-get clean

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
    && tar -xjf samtools-1.16.1.tar.bz2 \
    && rm samtools-1.16.1.tar.bz2 \
    && cd samtools-1.16.1 \
    && ./configure \
    && make \
    && make install \
    && mkdir -p /ROSE/src

RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

COPY src/* ./src