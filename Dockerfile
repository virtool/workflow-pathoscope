FROM rust:1.60.0-slim-buster as rust
WORKDIR /build
COPY /utils/eliminate_subtraction/ /build/
RUN cargo build -r

FROM virtool/workflow:4.2.2 as base
WORKDIR /app
COPY --from=rust /build/target/release/eliminate_subtraction ./
COPY workflow.py pathoscope.py ./

FROM virtool/workflow:4.2.2 as test
WORKDIR /test
COPY pyproject.toml poetry.lock ./
RUN pip install poetry
RUN poetry install
COPY --from=rust /build/target/release/eliminate_subtraction ./
COPY tests /test/tests
COPY workflow.py pathoscope.py ./
RUN ls
RUN poetry run pytest
