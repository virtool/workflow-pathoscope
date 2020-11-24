# Dockerfile for running `workflow-pathoscope` tests.

FROM workflow-pathoscope
WORKDIR test
COPY . .

RUN pip install tox

ENTRYPOINT ["tox"]

