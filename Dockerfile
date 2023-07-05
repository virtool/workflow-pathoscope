FROM python:3.10-buster as rust_utils
WORKDIR /build
RUN apt-get update && apt-get install -y curl build-essential
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN pip install maturin==0.14.5
COPY src src
COPY Cargo.toml Cargo.lock ./
RUN maturin build --release

FROM ghcr.io/virtool/workflow:5.3.1 as base
WORKDIR /app
RUN pip install --upgrade pip
COPY fixtures.py workflow.py pathoscope.py ./
COPY --from=rust_utils /build/target/wheels/rust_utils*.whl ./
RUN ls
RUN pip install rust_utils*.whl

FROM ghcr.io/virtool/workflow:5.3.1 as test
WORKDIR /test
RUN pip install --upgrade pip
COPY pyproject.toml poetry.lock ./
RUN curl -sSL https://install.python-poetry.org | python -
COPY tests /test/tests
COPY fixtures.py workflow.py pathoscope.py ./
COPY --from=rust_utils /build/target/wheels/rust_utils*.whl ./
RUN pip install rust_utils*.whl
RUN poetry install
RUN poetry add ./rust_utils*.whl
RUN ls
RUN poetry run pytest
