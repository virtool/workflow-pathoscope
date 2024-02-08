FROM debian:buster as prep
WORKDIR /build
RUN apt-get update && apt-get install -y make gcc zlib1g-dev wget unzip
RUN wget https://zlib.net/pigz/pigz-2.8.tar.gz && \
    tar -xzvf pigz-2.8.tar.gz && \
    cd pigz-2.8 && \
    make
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.2/bowtie2-2.3.2-legacy-linux-x86_64.zip && \
    unzip bowtie2-2.3.2-legacy-linux-x86_64.zip && \
    mkdir bowtie2 && \
    cp bowtie2-2.3.2-legacy/bowtie2* bowtie2


FROM python:3.10-buster as base
WORKDIR /app
COPY --from=prep /build/bowtie2/* /usr/local/bin/
COPY --from=prep /build/FastQC /opt/fastqc
COPY --from=prep /build/pigz-2.8/pigz /usr/local/bin/pigz
RUN chmod ugo+x /opt/fastqc/fastqc && \
    ln -fs /opt/fastqc/fastqc /usr/local/bin/fastqc && \
    for file in `ls /opt/hmmer/bin`; do ln -fs /opt/hmmer/bin/${file} /usr/local/bin/${file};  done
RUN apt-get update && \
    apt-get install -y --no-install-recommends curl build-essential default-jre && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -sSL https://install.python-poetry.org | python -
ENV PATH="/root/.cargo/bin:/root/.local/bin:${PATH}"
RUN pip install --upgrade pip
RUN pip install maturin==0.14.5
COPY src src
COPY Cargo.toml Cargo.lock fixtures.py utils.py  poetry.lock pyproject.toml workflow.py ./
RUN maturin build --release
RUN poetry export > requirements.txt
RUN pip install -r requirements.txt
RUN pip install /app/target/wheels/rust_utils*.whl
COPY VERSION* ./

FROM base as test
WORKDIR /app
RUN poetry export --with dev > requirements.txt
RUN pip install -r requirements.txt
COPY example ./example
COPY tests ./tests
