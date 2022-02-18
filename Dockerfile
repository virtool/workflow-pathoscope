FROM virtool/workflow:2.1.0 as base
WORKDIR /app
COPY workflow.py .
COPY pathoscope.py .
