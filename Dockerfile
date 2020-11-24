FROM virtool-workflow-standalone
COPY workflow.py .

ENTRYPOINT ["workflow", "--no-uv-loop", "run"]

