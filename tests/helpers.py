import os
import tempfile
import time
from pathlib import Path

from cichlid_data import Layout, LocalBackend, Store


class Env:
    """Two machines (a and b) sharing one fake cloud folder."""
    def __init__(self):
        self._tmp = tempfile.TemporaryDirectory()
        t = Path(self._tmp.name)
        self.cloud = t / "cloud"
        self.cloud.mkdir()
        self.a = Layout(t / "a", "fake:")
        self.b = Layout(t / "b", "fake:")

    def store(self, layout, **kw):
        return Store(layout, LocalBackend(self.cloud), **kw)

    def close(self):
        self._tmp.cleanup()


def write(path, text="x"):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


def touch_later(path):
    """Change a file's content and mtime so its fingerprint changes."""
    p = Path(path)
    p.write_text(p.read_text() + "+")
    t = time.time() + 5
    os.utime(p, (t, t))
