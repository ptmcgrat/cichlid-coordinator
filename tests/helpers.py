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


FAKE_RCLONE_PATH = str(Path(__file__).resolve().parent / "fake_rclone.py")


def coord_env(d):
    """Environment for code that uses the coordinator's settings: data root <d>/local, a fake cloud
    at <d>/cloud reached through the stand-in rclone, and no machine config file."""
    d = Path(d)
    return {"CICHLID_DATA_ROOT": str(d / "local"), "CICHLID_DATA_REMOTE": "fake:",
            "CICHLID_DATA_RCLONE": FAKE_RCLONE_PATH, "FAKE_CLOUD": str(d / "cloud"),
            "CICHLID_DATA_CONFIG": str(d / "no-config.toml")}


def releases_toml(parts, release="v1", status="building", deliver=("Mzeb",), inputs=None, description="test release"):
    """A releases.toml with one release. `parts` maps a part name to the text of its keys, e.g.
    {"TreeFam": 'queries = "project"', "OrthoFinder": 'species = ["Mzeb", "Drer"]'}."""
    text = f'[releases.{release}]\nstatus = "{status}"\ndescription = "{description}"\n'
    text += "deliver = [" + ", ".join(f'"{d}"' for d in deliver) + "]\n"
    text += f"[releases.{release}.inputs]\n" + "".join(f'{k} = "{v}"\n' for k, v in (inputs or {"treefam": "TreeFam 9"}).items())
    for part, body in parts.items():
        text += f"[releases.{release}.parts.{part}]\n{body}\n"
    return text
