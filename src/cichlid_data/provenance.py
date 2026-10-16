"""Who and what produced a transfer or a download; checksums."""
import hashlib
import os
import platform
import subprocess
import time

from . import __version__


def machine_name():
    """Readable machine name (on a Mac, the computer name: platform.node() is the VPN address on VPN)."""
    if platform.system() == "Darwin":
        try:
            r = subprocess.run(["scutil", "--get", "ComputerName"], capture_output=True, text=True, timeout=5)
            if r.returncode == 0 and r.stdout.strip():
                return r.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            pass
    return platform.node().split(".")[0] or platform.node()


def git_info(path):
    def git(*a):
        try:
            r = subprocess.run(["git", "-C", str(path)] + list(a), capture_output=True, text=True)
            return r.stdout.strip() if r.returncode == 0 else ""
        except OSError:
            return ""
    return {"branch": git("rev-parse", "--abbrev-ref", "HEAD") or "unknown",
            "commit": git("rev-parse", "--short", "HEAD") or "unknown",
            "uncommitted_changes": bool(git("status", "--porcelain", "--untracked-files=no"))}


def stamp(context="", repo=None):
    """Provenance of an upload: machine, user, time, coordinator version, and the caller's repo state."""
    out = {"context": context, "machine": machine_name(), "user": os.getenv("USER", ""),
           "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "coordinator": __version__}
    if repo:
        out.update(git_info(repo))
    return out


def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()
