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

FROM python:3.12.3-bookworm as deps
WORKDIR /app
COPY --from=bowtie2 /build/bowtie2/* /usr/local/bin/
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/fastqc /opt/fastqc
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /opt/hmmer /opt/hmmer
COPY --from=ghcr.io/virtool/workflow-tools:2.0.1 /usr/local/bin/pigz /usr/local/bin/
RUN apt-get update && \
    apt-get install -y --no-install-recommends default-jre && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

FROM python:3.12.3-bookworm as poetry
WORKDIR /app
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -sSL https://install.python-poetry.org | python -
ENV PATH="/root/.local/bin:/root/.cargo/bin:${PATH}" \
    POETRY_CACHE_DIR='/tmp/poetry_cache' \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1 \
    POETRY_VIRTUALENVS_CREATE=1
COPY Cargo.toml Cargo.lock poetry.lock pyproject.toml ./
COPY python/ ./python
COPY src ./src
RUN poetry install
RUN poetry run maturin develop --release

FROM deps as base
WORKDIR /app
ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:/opt/fastqc:/opt/hmmer/bin:${PATH}"
RUN chmod ugo+x /opt/fastqc/fastqc
COPY --from=poetry /app/.venv /app/.venv
COPY --from=poetry /app/python /app/python
COPY fixtures.py workflow.py VERSION* ./

FROM deps as test
WORKDIR /app
ENV PATH="/root/.local/bin:/root/.cargo/bin:${PATH}"
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -sSL https://install.python-poetry.org | python -
COPY Cargo.lock Cargo.toml pyproject.toml poetry.lock ./
COPY src ./src
COPY python ./python
RUN poetry install
RUN poetry run maturin develop --release
COPY example ./example
COPY tests ./tests
COPY fixtures.py workflow.py VERSION* ./
