FROM python:3.12.3-bookworm AS deps
WORKDIR /app
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/bowtie2/2.5.4/bowtie* /usr/local/bin/
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/hmmer/3.2.1 /opt/hmmer
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/fastqc/0.11.9 /opt/fastqc
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/pigz/2.8/pigz /usr/local/bin/
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/samtools/1.22.1/bin/samtools /usr/local/bin/
RUN apt-get update && \
    apt-get install -y --no-install-recommends default-jre && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean

FROM python:3.12.3-bookworm AS uv
WORKDIR /app
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
RUN apt-get update && \
    apt-get install -y --no-install-recommends libclang-dev && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean
ENV PATH="/root/.cargo/bin:/root/.local/bin:${PATH}" \
    UV_CACHE_DIR='/tmp/uv_cache'
COPY Cargo.toml Cargo.lock uv.lock pyproject.toml README.md ./
COPY python/ ./python
COPY src ./src
RUN uv sync
RUN uv run maturin develop --release

FROM deps AS base
WORKDIR /app
ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:/opt/fastqc:/opt/hmmer/bin:${PATH}"
RUN chmod ugo+x /opt/fastqc/fastqc
COPY --from=uv /app/.venv /app/.venv
COPY --from=uv /app/python /app/python
COPY fixtures.py workflow.py VERSION* ./

FROM deps AS dev
WORKDIR /app
ENV PATH="/root/.cargo/bin:/root/.local/bin:${PATH}" \
    UV_CACHE_DIR='/tmp/uv_cache'
RUN apt-get update && \
    apt-get install -y --no-install-recommends libclang-dev && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get clean
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
RUN curl -LsSf https://astral.sh/uv/install.sh | sh

FROM dev AS test
COPY Cargo.lock Cargo.toml pyproject.toml uv.lock ./
COPY README.md ./
COPY src ./src
COPY python ./python
RUN uv sync
RUN uv run maturin develop --release
COPY example ./example
COPY tests ./tests
COPY fixtures.py workflow.py VERSION* ./
