import subprocess
import pytest

def test_help_runs():
    """Ensure that --help executes without error."""
    result = subprocess.run(["python", "-m", "structgen.cli", "--help"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "usage:" in result.stdout.lower()
