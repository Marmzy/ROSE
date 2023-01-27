FROM python:3.8.10

WORKDIR /ROSE

COPY ROSE.sh requirements.txt README.md LICENSE ./

RUN apt-get update && \
    apt-get install -y \
    gcc \
    make \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \ 
    libncursesw5-dev \
    liblzma-dev

RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
    tar -xjf samtools-1.16.1.tar.bz2 && \
    rm samtools-1.16.1.tar.bz2 && \
    cd samtools-1.16.1 && \
    ./configure && \
    make && \
    make install && \
    mkdir -p /ROSE/src

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

COPY src/* ./src