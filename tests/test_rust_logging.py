"""Test Rust logging integration with Python logging systems."""

import os

import pytest

from workflow_pathoscope.rust import init_logging


def test_init_logging_basic():
    """Test that init_logging can be called without errors."""
    init_logging("info")
    init_logging("debug")
    init_logging("warn")
    init_logging("error")
    init_logging("trace")


def test_init_logging_with_none():
    """Test init_logging with None parameter uses default level."""
    init_logging(None)


@pytest.mark.parametrize("level", ["trace", "debug", "info", "warn", "error"])
def test_logging_levels(level):
    """Test different logging levels."""
    init_logging(level)


def test_multiple_init_calls():
    """Test that multiple calls to init_logging are safe."""
    init_logging("info")
    init_logging("debug")
    init_logging("warn")


def test_env_var_parsing():
    """Test RUST_LOG environment variable support."""

    original_rust_log = os.environ.get("RUST_LOG")

    try:
        os.environ["RUST_LOG"] = "debug"
        init_logging(None)

        if "RUST_LOG" in os.environ:
            del os.environ["RUST_LOG"]

        init_logging(None)

    finally:
        if original_rust_log is not None:
            os.environ["RUST_LOG"] = original_rust_log

        elif "RUST_LOG" in os.environ:
            del os.environ["RUST_LOG"]


def test_logging_output_capture(log):
    """Test that Rust logs can be captured by Python logging."""
    init_logging("info")

    assert log.has("rust logging initialized")
