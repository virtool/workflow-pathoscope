FROM virtool/workflow:2.1.3 as base
WORKDIR /app
COPY workflow.py .
COPY pathoscope.py .
