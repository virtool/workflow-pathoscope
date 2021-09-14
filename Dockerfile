FROM virtool/workflow:nightly as base
WORKDIR /app
COPY workflow.py .
COPY pathoscope.py .

FROM base as test
WORKDIR /app
RUN ["pip", "install", "poetry"]
COPY poetry.lock .
COPY pyproject.toml .
RUN ["poetry", "install"]
COPY tests tests
RUN ["ls", "tests"]
RUN ["poetry", "run", "pytest", "-x"]
