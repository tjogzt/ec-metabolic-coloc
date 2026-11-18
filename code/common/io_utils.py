"""Input-path and remote-file helpers shared by colocalisation stages."""

from __future__ import annotations

import pathlib
import subprocess

from config.config import PATHS


def resolve_input(name: str) -> pathlib.Path:
    """Resolve a user-supplied file name against the data directory."""
    path = pathlib.Path(name)
    return path if path.is_absolute() else PATHS.data / path


def remote_size(url: str) -> int:
    """Content-Length of a remote file, or 0 when it cannot be read.

    Zero is not treated as an error: the run may legitimately be offline, in
    which case the local file cannot be validated but can still be used.
    """
    proc = subprocess.run(["curl", "-sI", "-m", "60", url],
                          capture_output=True, text=True)
    for line in proc.stdout.splitlines():
        if line.lower().startswith("content-length"):
            try:
                return int(line.split(":", 1)[1].strip())
            except ValueError:
                return 0
    return 0
