FROM virtool/workflow:1.0.1 as base
WORKDIR /app
COPY workflow.py .
COPY pathoscope.py .
