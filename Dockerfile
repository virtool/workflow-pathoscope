FROM virtool/workflow:2.1.2 as base
WORKDIR /app
COPY workflow.py .
COPY pathoscope.py .
