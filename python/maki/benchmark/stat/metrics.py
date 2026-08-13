from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from tempfile import NamedTemporaryFile

from .schema import BenchmarkMetrics


class CommandExecutionError(RuntimeError):
    def __init__(
        self,
        command: list[str],
        returncode: int,
        stdout: str,
        stderr: str,
    ):
        super().__init__(
            f"Command failed with exit code {returncode}: "
            f"{' '.join(command)}\n"
            f"STDERR:\n{stderr}"
        )
        self.command = command
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


class TimeWrapper:
    """
    Wrapper around GNU /usr/bin/time.
    """

    _FORMAT = json.dumps(
        {
            "elapsed_seconds": "%e",
            "user_seconds": "%U",
            "system_seconds": "%S",
            "cpu_percent": "%P",
            "max_rss_kb": "%M",
            "major_page_faults": "%F",
            "minor_page_faults": "%R",
            "voluntary_context_switches": "%w",
            "involuntary_context_switches": "%c",
            "filesystem_inputs": "%I",
            "filesystem_outputs": "%O",
            "exit_code": "%x",
        }
    )

    def __init__(self, time_path: str = "/usr/bin/time"):
        self.time_path = Path(time_path)

        if not self.time_path.exists():
            raise FileNotFoundError(
                f"{self.time_path} does not exist"
            )

        if not self.time_path.is_file():
            raise ValueError(
                f"{self.time_path} is not a regular file"
            )

        if not shutil.which(str(self.time_path)):
            raise RuntimeError(
                f"{self.time_path} is not executable"
            )

    def run(
        self,
        command: list[str],
        *,
        raise_on_failure: bool = True
    ) -> BenchmarkMetrics:
        """
        Run a command and collect timing statistics.

        Example:
            stats = wrapper.run(["sleep", "2"])
        """

        with NamedTemporaryFile(mode="r+") as stats_file:
            try:
                proc = subprocess.run(
                    [
                        str(self.time_path),
                        "-f",
                        self._FORMAT,
                        "-o",
                        stats_file.name,
                        *command,
                    ],
                    text=True,
                    capture_output=True,
                    check=False
                )

            except FileNotFoundError as exc:
                raise RuntimeError(
                    f"Unable*to execute command: {command[0]!r}*"
                    f"was not found") from exc

            except PermissionError as exc:
              raise RuntimeError(
                   f"Permission denied executin* command: "
                    f"{command[0]!r}"
                ) from exc
              
            except OSError as exc:
                 raise RuntimeError(
                      f"Failed to execute comma*d "
                        f"{' '.join(command)}: {exc}"
                    ) from exc

            stats_file.seek(0)
            raw = stats_file.read().strip()

        json_line = None
        
        for line in reversed(raw.splitlines()):
            line = line.strip()
            if line.startswith("{") and line.endswith("}"):
                json_line = line
                break
            
        if json_line is None:
            raise RuntimeError(
                f"Failed to locate JSON output from {self.time_path}\n"
                f"Raw output was:\n{raw}"
            )
        
        data = json.loads(json_line)

        stats = BenchmarkMetrics(
            elapsed_seconds=float(data["elapsed_seconds"]),
            user_seconds=float(data["user_seconds"]),
            system_seconds=float(data["system_seconds"]),
            cpu_percent=float(data["cpu_percent"].rstrip("%")),
            max_rss_kb=int(data["max_rss_kb"]),
            major_page_faults=int(data["major_page_faults"]),
            minor_page_faults=int(data["minor_page_faults"]),
            voluntary_context_switches=int(
                data["voluntary_context_switches"]
            ),
            involuntary_context_switches=int(
                data["involuntary_context_switches"]
            ),
            filesystem_inputs=int(data["filesystem_inputs"]),
            filesystem_outputs=int(data["filesystem_outputs"]),
            exit_code=int(data["exit_code"]),
        )
        
        if raise_on_failure and proc.returncode != 0:
            raise CommandExecutionError(
                command=command,
                returncode=proc.returncode,
                stdout=proc.stdout,
                stderr=proc.stderr,
            )
        
        return stats
        