"""Test Rust logging integration with Python logging systems."""

import os

import pytest

from workflow_pathoscope.rust import init_logging


def test_init_logging_basic():
    """Test that init_logging can be called without errors."""
    # This should not raise any exceptions
    init_logging("info")
    init_logging("debug")
    init_logging("warn")
    init_logging("error")
    init_logging("trace")


def test_init_logging_with_none():
    """Test init_logging with None parameter."""
    # Should use default level
    init_logging(None)


@pytest.mark.parametrize("level", ["trace", "debug", "info", "warn", "error"])
def test_logging_levels(level):
    """Test different logging levels."""
    init_logging(level)


def test_multiple_init_calls():
    """Test that multiple calls to init_logging are safe."""
    # Should be safe to call multiple times due to std::sync::Once
    init_logging("info")
    init_logging("debug")
    init_logging("warn")


def test_env_var_parsing():
    """Test RUST_LOG environment variable support."""

    # Save original value
    original_rust_log = os.environ.get("RUST_LOG")

    try:
        # Test with RUST_LOG set
        os.environ["RUST_LOG"] = "debug"
        init_logging(None)  # Should use env var

        # Test with RUST_LOG unset
        if "RUST_LOG" in os.environ:
            del os.environ["RUST_LOG"]
        init_logging(None)  # Should use default

    finally:
        # Restore original value
        if original_rust_log is not None:
            os.environ["RUST_LOG"] = original_rust_log
        elif "RUST_LOG" in os.environ:
            del os.environ["RUST_LOG"]


def test_logging_output_capture(log):
    """Test that Rust logs can be captured by Python logging."""
    # Initialize Rust logging
    init_logging("info")

    # Check that the Rust logging initialization was captured
    assert log.has("Rust logging initialized with custom Python bridge")


def test_structlog_detection(log):
    """Test that structlog is properly detected when available."""
    init_logging("info")

    # Verify structlog captured the log output
    assert log.has("Rust logging initialized with custom Python bridge")
