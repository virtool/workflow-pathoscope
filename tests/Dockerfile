FROM virtool/external-tools:0.2.0

WORKDIR test

RUN pip install poetry
COPY pyproject.toml .
COPY poetry.lock .
RUN poetry install

COPY . .

ENTRYPOINT ["poetry", "run", "pytest"]
