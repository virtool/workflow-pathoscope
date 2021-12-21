FROM virtool/workflow:1.1.0 as base
WORKDIR /app
COPY workflow.py .
COPY pathoscope.py .
