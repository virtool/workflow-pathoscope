FROM debian:bookworm as bowtie2
WORKDIR /build
RUN apt-get update && apt-get install -y build-essential cmake wget zlib1g-dev
RUN wget https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.5.4.tar.gz
RUN tar -xvf v2.5.4.tar.gz
WORKDIR bowtie2-2.5.4
RUN make
RUN mkdir /build/bowtie2
RUN cp bowtie2* /build/bowtie2/

FROM debian:bookworm as pigz
WORKDIR /build
RUN apt-get update && apt-get install -y gcc make wget zlib1g-dev
RUN wget https://zlib.net/pigz/pigz-2.8.tar.gz && \
    tar -xzvf pigz-2.8.tar.gz && \
    cd pigz-2.8 && \
    make

FROM python:3.12.3-bookworm
COPY --from=bowtie2 /build/bowtie2/* /usr/local/bin/
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/fastqc /opt/fastqc
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/hmmer /opt/hmmer
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/pigz /usr/local/bin/
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -sSL https://install.python-poetry.org | python -
ENV PATH="/root/.local/bin:/root/.cargo/bin:${PATH}"
