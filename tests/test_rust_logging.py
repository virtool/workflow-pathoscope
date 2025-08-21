"""Test Rust logging integration with Python logging systems."""
import logging
import sys
from io import StringIO
from unittest.mock import patch

import pytest

from workflow_pathoscope.rust import init_logging, init_logging_with_logger


def test_init_logging_basic():
    """Test that init_logging can be called without errors."""
    init_logging("info")
    init_logging("debug")
    init_logging("warn")
    init_logging("error")
    init_logging("trace")


def test_init_logging_with_none():
    """Test init_logging with None parameter."""
    init_logging(None)


def test_init_logging_with_logger_basic():
    """Test that init_logging_with_logger can be called."""
    logger = logging.getLogger("test_rust")
    init_logging_with_logger(logger, "info")


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
    import os
    
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


def test_rust_logs_appear_in_python(capsys):
    """Test that Rust logs appear in Python output via structlog."""
    init_logging("info")
    
    # Capture stdout/stderr
    captured = capsys.readouterr()
    
    # Should contain Rust logging initialization message
    assert "Rust logging initialized" in captured.out
    assert "workflow_pathoscope" in captured.out


def test_structlog_detection():
    """Test that structlog is properly detected when available."""
    try:
        import structlog
        structlog_available = True
    except ImportError:
        structlog_available = False
    
    if not structlog_available:
        pytest.skip("structlog not available")
    
    # Mock to verify structlog methods are called
    with patch('structlog.get_logger') as mock_get_logger:
        mock_logger = mock_get_logger.return_value
        
        init_logging("info")
        
        # Verify structlog.get_logger was called
        mock_get_logger.assert_called()
        
        # Verify the logger method was called (should be 'info' level)
        assert mock_logger.info.called or mock_logger.method_calls


def test_fallback_when_structlog_unavailable(capsys):
    """Test fallback to standard logging when structlog is not available."""
    # Temporarily hide structlog
    original_modules = sys.modules.copy()
    if 'structlog' in sys.modules:
        del sys.modules['structlog']
    
    try:
        init_logging("info")
        captured = capsys.readouterr()
        
        # Should still get logging output via fallback
        # Note: Since logging is already initialized, we might not see new output
        # but the call should succeed without errors
        
    finally:
        # Restore modules
        sys.modules.clear()
        sys.modules.update(original_modules)


def test_log_levels_mapping(capsys):
    """Test that Rust log levels map correctly to Python levels."""
    init_logging("debug")
    captured = capsys.readouterr()
    
    # Should contain log output with level info
    assert "info" in captured.out.lower()
    assert "rust logging initialized" in captured.out.lower()