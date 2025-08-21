"""
Test logging integration between Rust and Python using our custom logging bridge.
"""
import logging
import sys
from io import StringIO

import pytest

from workflow_pathoscope.rust import init_logging, init_logging_with_logger


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


def test_init_logging_with_logger_basic():
    """Test that init_logging_with_logger can be called."""
    logger = logging.getLogger("test_rust")
    # This should not raise any exceptions
    init_logging_with_logger(logger, "info")


def test_logging_output_capture():
    """Test that Rust logs can be captured by Python logging."""
    # Set up a string buffer to capture log output
    log_capture = StringIO()
    handler = logging.StreamHandler(log_capture)
    formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    
    # Configure Python logging
    root_logger = logging.getLogger()
    original_level = root_logger.level
    original_handlers = root_logger.handlers[:]
    
    # Clear existing handlers and add our test handler
    root_logger.handlers.clear()
    root_logger.addHandler(handler)
    root_logger.setLevel(logging.INFO)
    
    try:
        # Initialize Rust logging
        init_logging("info")
        
        # Get the captured output
        output = log_capture.getvalue()
        
        # Should contain the initialization message from Rust
        assert "Rust logging initialized with custom Python bridge" in output
        assert "workflow_pathoscope.logging" in output or "logging" in output
        
    finally:
        # Clean up - restore original logging configuration
        root_logger.handlers.clear()
        for handler in original_handlers:
            root_logger.addHandler(handler)
        root_logger.setLevel(original_level)


@pytest.mark.parametrize("level", ["trace", "debug", "info", "warn", "error"])
def test_logging_levels(level):
    """Test different logging levels."""
    init_logging(level)
    # Should not raise any exceptions
    assert True


def test_multiple_init_calls():
    """Test that multiple calls to init_logging are safe."""
    # Should be safe to call multiple times due to std::sync::Once
    init_logging("info")
    init_logging("debug")
    init_logging("warn")
    # Should not cause any issues
    assert True


def test_rust_log_levels_mapping():
    """Test that different log levels work as expected."""
    # Set up capture for all levels
    log_capture = StringIO()
    handler = logging.StreamHandler(log_capture)
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    handler.setFormatter(formatter)
    
    # Configure Python logging to capture all levels
    root_logger = logging.getLogger()
    original_level = root_logger.level
    original_handlers = root_logger.handlers[:]
    
    root_logger.handlers.clear()
    root_logger.addHandler(handler)
    root_logger.setLevel(logging.DEBUG)  # Capture all levels
    
    try:
        # Initialize Rust logging with debug level
        init_logging("debug")
        
        # Get the output - should contain the initialization message
        output = log_capture.getvalue()
        assert "INFO: Rust logging initialized" in output
        
    finally:
        # Clean up
        root_logger.handlers.clear()
        for handler in original_handlers:
            root_logger.addHandler(handler)
        root_logger.setLevel(original_level)


def test_env_var_parsing():
    """Test RUST_LOG environment variable support."""
    import os
    
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


def test_structlog_integration():
    """Test that Rust logging integrates properly with structlog."""
    try:
        import structlog
    except ImportError:
        pytest.skip("structlog not available")
    
    from io import StringIO
    import sys
    
    # Capture output
    output_buffer = StringIO()
    
    # Configure structlog with simple console output
    structlog.configure(
        processors=[
            structlog.stdlib.add_log_level,
            structlog.stdlib.add_logger_name,
            structlog.dev.ConsoleRenderer(colors=False)
        ],
        wrapper_class=structlog.stdlib.BoundLogger,
        logger_factory=structlog.WriteLoggerFactory(file=output_buffer),
        context_class=dict,
        cache_logger_on_first_use=True,
    )
    
    # Initialize Rust logging
    init_logging("info")
    
    # Get the output
    output = output_buffer.getvalue()
    
    # Should contain structured log output from Rust initialization
    # The exact format depends on structlog configuration but should be more
    # than just plain text
    assert len(output) > 0
    print(f"Structlog output: {output}")  # For debugging


def test_fallback_to_standard_logging():
    """Test fallback when structlog is not available."""
    # Temporarily hide structlog module
    import sys
    structlog_module = sys.modules.get('structlog')
    if 'structlog' in sys.modules:
        del sys.modules['structlog']
    
    try:
        # Set up capture for standard logging
        log_capture = StringIO()
        handler = logging.StreamHandler(log_capture)
        formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        
        # Configure Python logging
        root_logger = logging.getLogger()
        original_level = root_logger.level
        original_handlers = root_logger.handlers[:]
        
        # Clear existing handlers and add our test handler
        root_logger.handlers.clear()
        root_logger.addHandler(handler)
        root_logger.setLevel(logging.INFO)
        
        try:
            # Initialize Rust logging (should fall back to standard logging)
            init_logging("info")
            
            # Get the captured output
            output = log_capture.getvalue()
            
            # Should contain the initialization message
            assert "Rust logging initialized" in output
            
        finally:
            # Clean up - restore original logging configuration
            root_logger.handlers.clear()
            for h in original_handlers:
                root_logger.addHandler(h)
            root_logger.setLevel(original_level)
    
    finally:
        # Restore structlog module if it was available
        if structlog_module is not None:
            sys.modules['structlog'] = structlog_module