# Build & Test Instructions

## Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_specific.py -v

# Run with coverage
pytest tests/ --cov=src --cov-report=term
```

## Environment

This project requires Python 3.9+ with dependencies in `pyproject.toml`.

Modal is used for GPU compute but not required for running tests locally.

## Git Workflow

After making fixes:
```bash
git add <modified_files>
git commit -m "fix: descriptive message"
```
