//! Custom PyO3 logging bridge implementation
//! 
//! This module provides a custom logging bridge that forwards Rust log messages
//! to Python's logging system. Inspired by pyo3-log but simplified and tailored
//! for our specific needs.

use pyo3::prelude::*;
use pyo3::types::{PyTuple, PyDict};
use log::{Level, LevelFilter, Log, Metadata, Record};
use std::sync::Once;

/// Environment variable name for Rust log level configuration
static RUST_LOG: &str = "RUST_LOG";

/// Custom logger that forwards Rust logs to Python
#[derive(Debug)]
struct PythonLogger {
    filter: LevelFilter,
}

impl PythonLogger {
    fn new(filter: LevelFilter) -> Self {
        Self { filter }
    }
}

impl Log for PythonLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= self.filter
    }

    fn log(&self, record: &Record) {
        if !self.enabled(record.metadata()) {
            return;
        }

        // Forward to Python with GIL management
        let _ = Python::with_gil(|py| {
            forward_to_python(py, record)
        });
    }

    fn flush(&self) {}
}

/// Forward a Rust log record to Python's logging system
/// 
/// This function attempts to use structlog first for structured logging,
/// then falls back to standard Python logging if structlog is not available.
fn forward_to_python(py: Python, record: &Record) -> PyResult<()> {
    let target = "rust";

    // Format the log message
    let message = format!("{}", record.args());
    
    // Try to use structlog first for structured logging
    if let Ok(structlog) = py.import_bound("structlog") {
        if let Ok(get_logger) = structlog.getattr("get_logger") {
            // Get structlog bound logger
            let logger = get_logger.call1((target,))?;
            
            // Create structured data as keyword arguments
            let kwargs = PyDict::new_bound(py);
            kwargs.set_item("event", &message)?;
            
            // Call appropriate structlog method with structured data
            let method_name = match record.level() {
                Level::Error => "error",
                Level::Warn => "warning", 
                Level::Info => "info",
                Level::Debug => "debug",
                Level::Trace => "debug", // structlog doesn't have trace, use debug
            };
            
            // Call the structlog method with kwargs
            logger.call_method(method_name, (), Some(&kwargs))?;
            return Ok(());
        }
    }
    
    // Fallback to standard Python logging if structlog is not available
    fallback_to_standard_logging(py, record, &target, &message)
}

/// Fallback function for standard Python logging when structlog is not available
fn fallback_to_standard_logging(
    py: Python, 
    record: &Record, 
    target: &str, 
    message: &str
) -> PyResult<()> {
    // Import Python logging module
    let logging = py.import_bound("logging")?;
    
    // Get Python logger for this target
    let logger = logging.call_method1("getLogger", (target,))?;
    
    // Map Rust log level to Python numeric level
    let py_level = rust_level_to_python(record.level());
    
    // Check if logger would log at this level (for performance)
    let is_enabled: bool = logger.call_method1("isEnabledFor", (py_level,))?.extract()?;
    
    if is_enabled {
        // Create a proper log record using Python's makeRecord method
        let none = py.None();
        let log_record = logger.call_method1(
            "makeRecord",
            (
                target,                          // name
                py_level,                        // level
                record.file().unwrap_or(""),     // pathname
                record.line().unwrap_or(0),      // lineno
                message,                         // msg
                PyTuple::empty_bound(py),        // args
                &none,                           // exc_info
            ),
        )?;
        
        // Send the record to Python logging handlers
        logger.call_method1("handle", (&log_record,))?;
    }
    
    Ok(())
}

/// Map Rust log levels to Python numeric levels
fn rust_level_to_python(level: Level) -> u8 {
    match level {
        Level::Error => 40,  // logging.ERROR
        Level::Warn => 30,   // logging.WARNING  
        Level::Info => 20,   // logging.INFO
        Level::Debug => 10,  // logging.DEBUG
        Level::Trace => 5,   // Custom level between DEBUG and NOTSET
    }
}

/// Initialize Rust logging to forward to Python's logging system
/// 
/// This function should be called from Python to set up logging for Rust components.
/// It installs a custom logger that forwards all Rust log messages to Python's
/// logging system.
/// 
/// # Arguments
/// * `log_level` - Optional log level string ("debug", "info", "warn", "error", "trace")
///                If None, defaults to "info" or uses RUST_LOG environment variable
/// 
/// # Examples
/// ```python
/// import workflow_pathoscope.rust as rust
/// import logging
/// logging.basicConfig(level=logging.INFO)
/// rust.init_logging("info")
/// ```
#[pyfunction]
#[pyo3(signature = (log_level=None))]
pub fn init_logging(_py: Python, log_level: Option<String>) -> PyResult<()> {
    static INIT: Once = Once::new();
    
    INIT.call_once(|| {
        // Determine log level filter
        let level_filter = if let Some(level) = log_level {
            match level.to_lowercase().as_str() {
                "trace" => LevelFilter::Trace,
                "debug" => LevelFilter::Debug,
                "info" => LevelFilter::Info,
                "warn" => LevelFilter::Warn,
                "error" => LevelFilter::Error,
                _ => LevelFilter::Info,
            }
        } else {
            // Use RUST_LOG environment variable or default to info
            if let Ok(rust_log) = std::env::var(RUST_LOG) {
                // Simple parsing for common cases
                match rust_log.to_lowercase().as_str() {
                    "trace" => LevelFilter::Trace,
                    "debug" => LevelFilter::Debug,
                    "info" => LevelFilter::Info,
                    "warn" => LevelFilter::Warn,
                    "error" => LevelFilter::Error,
                    _ => LevelFilter::Info,
                }
            } else {
                LevelFilter::Info
            }
        };
        
        // Create and install our custom logger
        let logger = PythonLogger::new(level_filter);
        
        if log::set_boxed_logger(Box::new(logger)).is_err() {
            // Logger already installed, ignore
            eprintln!("Warning: Logger already installed, ignoring init_logging call");
        }
        
        // Set maximum log level for performance
        log::set_max_level(level_filter);
    });
    
    // Test that logging works
    log::info!("Rust logging initialized with custom Python bridge");
    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_rust_level_to_python() {
        assert_eq!(rust_level_to_python(Level::Error), 40);
        assert_eq!(rust_level_to_python(Level::Warn), 30);
        assert_eq!(rust_level_to_python(Level::Info), 20);
        assert_eq!(rust_level_to_python(Level::Debug), 10);
        assert_eq!(rust_level_to_python(Level::Trace), 5);
    }
    
    #[test]
    fn test_target_conversion() {
        let target = "workflow_pathoscope::matrix::build_matrix";
        let expected = "workflow_pathoscope.matrix.build_matrix";
        assert_eq!(target.replace("::", "."), expected);
    }
    
    #[test]
    fn test_init_logging_basic() {
        // Prepare Python for testing
        pyo3::prepare_freethreaded_python();
        
        Python::with_gil(|py| {
            let result = init_logging(py, Some("info".to_string()));
            assert!(result.is_ok());
        });
    }
}
