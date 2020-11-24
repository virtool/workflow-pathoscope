# Dockerfile for the `pathoscope` Virtool workflow.

# Bowtie2
FROM alpine:latest as bowtie
WORKDIR /build
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.2/bowtie2-2.3.2-legacy-linux-x86_64.zip && \
    unzip bowtie2-2.3.2-legacy-linux-x86_64.zip && \
    mkdir bowtie2 && \
    cp bowtie2-2.3.2-legacy/bowtie2* bowtie2


FROM virtool-workflow-standalone
COPY --from=bowtie /build/bowtie2/* /usr/local/bin/
COPY workflow.py .

ENTRYPOINT ["workflow", "--no-uv-loop", "run"]

