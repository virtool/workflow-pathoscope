FROM virtool/workflow:1.0.8 as base
WORKDIR /app
COPY workflow.py .
COPY pathoscope.py .
