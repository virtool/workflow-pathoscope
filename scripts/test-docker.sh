#!/bin/bash
set -euo pipefail

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

IMAGE_NAME="pathoscope:dev"
CONTAINER_NAME="pathoscope-test-$$"

# Function to cleanup on exit
cleanup() {
    if docker ps -q -f name="$CONTAINER_NAME" | grep -q .; then
        echo -e "${YELLOW}Cleaning up container...${NC}"
        docker rm -f "$CONTAINER_NAME" >/dev/null 2>&1 || true
    fi
}

# Set trap to cleanup on script exit
trap cleanup EXIT INT TERM

# Check if image exists, build if not
if ! docker image inspect "$IMAGE_NAME" >/dev/null 2>&1; then
    echo -e "${BLUE}Building dev image...${NC}"
    docker build --target dev -t "$IMAGE_NAME" .
fi

# Parse command line arguments
MODE="test"
PYTEST_ARGS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --watch)
            MODE="watch"
            shift
            ;;
        --shell)
            MODE="shell"
            shift
            ;;
        *)
            PYTEST_ARGS="$PYTEST_ARGS $1"
            shift
            ;;
    esac
done

case $MODE in
    "shell")
        echo -e "${GREEN}Starting interactive shell...${NC}"
        docker run -it --rm \
            --name "$CONTAINER_NAME" \
            -v "$(pwd):/app" \
            -w /app \
            "$IMAGE_NAME" \
            bash
        ;;
    "watch")
        echo -e "${GREEN}Starting test watch mode...${NC}"
        docker run -it --rm \
            --name "$CONTAINER_NAME" \
            -v "$(pwd):/app" \
            -w /app \
            "$IMAGE_NAME" \
            bash -c "poetry install --no-root && poetry run maturin develop --release && poetry run pytest-watch$PYTEST_ARGS"
        ;;
    "test")
        echo -e "${GREEN}Running tests...${NC}"
        docker run --rm \
            --name "$CONTAINER_NAME" \
            -v "$(pwd):/app" \
            -w /app \
            "$IMAGE_NAME" \
            bash -c "echo 'Running Python tests...' && poetry install --no-root && poetry run maturin develop --release && poetry run pytest$PYTEST_ARGS"
        ;;
esac