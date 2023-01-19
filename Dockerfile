FROM rust:1.60.0-slim-buster as rust
WORKDIR /build
COPY /utils/eliminate_subtraction/ /build/
RUN cargo build -r

FROM python:3.10-buster as rustExpectMax
WORKDIR /build
RUN apt-get update && apt-get install -y curl build-essential
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN pip install maturin==0.14.5
COPY src src
COPY Cargo.toml Cargo.lock ./
RUN maturin build --release

FROM virtool/workflow:5.2.1 as base
WORKDIR /app
RUN pip install --upgrade pip
COPY pyproject.toml poetry.lock ./
COPY --from=rust /build/target/release/eliminate_subtraction ./
COPY fixtures.py workflow.py pathoscope.py ./
COPY --from=rustExpectMax /build/target/wheels/virtool_expectation_maximization*.whl ./
RUN pip install virtool_expectation_maximization*.whl
RUN poetry install
RUN poetry add ./virtool_expectation_maximization*.whl

FROM virtool/workflow:5.2.1 as test
WORKDIR /test
RUN pip install --upgrade pip
COPY pyproject.toml poetry.lock ./
RUN curl -sSL https://install.python-poetry.org | python -
COPY --from=rust /build/target/release/eliminate_subtraction ./
COPY tests /test/tests
COPY fixtures.py workflow.py pathoscope.py ./
COPY --from=rustExpectMax /build/target/wheels/virtool_expectation_maximization*.whl ./
RUN pip install virtool_expectation_maximization*.whl
RUN poetry install
RUN poetry add ./virtool_expectation_maximization*.whl
RUN ls
RUN poetry run pytest