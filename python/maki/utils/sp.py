import subprocess
import logging
import sys

logger = logging.getLogger(__name__)


def run_in_current_env(cmd):
    cmd = [
        sys.executable, "-m",
        *cmd
    ]
    
    try:
        return subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        logger.error("Command failed: %s", " ".join(map(str, e.cmd)))
        logger.error("stderr:\n%s", e.stderr)
        raise
