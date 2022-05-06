FROM virtool/workflow:4.0.2 as base
WORKDIR /app
COPY workflow.py .
COPY pathoscope.py .
