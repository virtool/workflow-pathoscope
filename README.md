# workflow-pathoscope

A workflow for detecting viruses in high-throughput sequencing (HTS) libraries.

## Contributing

### Run Tests Using docker-compose

The tests require MongoDB and Redis to be available in order to run. A 
[tests/docker-compose.yml](tests/docker-compose.yml) file is available to
run the test suite alongside these services. 

To run the tests in docker-compose use;

```shell script
docker build -t workflow-pathoscope .
cd tests
docker-compose build
docker-compose up --abort-on-container-exit 
```


You may need to prefix these commands with `sudo` depending on your environment.

### Run The Workflow Alongside Mongo and Redis

Another [docker-compose.yml](docker-compose.yml) is available for running the workflow
container with Mongo and Redis. 

Before running docker-compose, ensure the `virtool-workflow-standalone` image is
available. 

```shell script
docker build -t virtool-workflow-standalone https://raw.githubusercontent.com/virtool/virtool-workflow/master/docker/Dockerfile
```

From the repository root;

```shell script
docker-compose build 
docker-compose up
```
