"""Settings: where this machine keeps the data, and where the cloud copy lives.

Read from ~/.config/cichlid-data/config.toml (or the file named by CICHLID_DATA_CONFIG), with
environment variables taking precedence:

    local_root   = "~/Temp/TreeFamData"                                   # CICHLID_DATA_ROOT
    cloud_remote = "ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/TreeFamData"  # CICHLID_DATA_REMOTE
    rclone       = "rclone"                                               # CICHLID_DATA_RCLONE

`cichlid-data init` writes the file with these defaults (on utaka, local_root defaults to
/Data/$USER/Temp/TreeFamData).
"""
import os
import platform
import tomllib
from dataclasses import dataclass
from pathlib import Path

DEFAULT_REMOTE = "ptm_dropbox:/COS/BioSci/BioSci-McGrath/Apps/TreeFamData"


class ConfigError(RuntimeError):
    pass


def config_dir():
    """The coordinator's config/ folder (species list, species sets, sources): the repository's own,
    unless CICHLID_DATA_CONFIG_DIR points elsewhere."""
    return Path(os.environ.get("CICHLID_DATA_CONFIG_DIR") or Path(__file__).resolve().parents[2] / "config")


def config_path():
    return Path(os.environ.get("CICHLID_DATA_CONFIG", "~/.config/cichlid-data/config.toml")).expanduser()


def default_local_root():
    if "utaka" in platform.node():
        return f"/Data/{os.getenv('USER', 'user')}/Temp/TreeFamData"
    return "~/Temp/TreeFamData"


@dataclass(frozen=True)
class Settings:
    local_root: Path
    cloud_remote: str
    rclone: str = "rclone"

    @classmethod
    def load(cls):
        values = {}
        p = config_path()
        if p.is_file():
            values = tomllib.loads(p.read_text())
        for key, env in (("local_root", "CICHLID_DATA_ROOT"), ("cloud_remote", "CICHLID_DATA_REMOTE"),
                         ("rclone", "CICHLID_DATA_RCLONE")):
            if os.environ.get(env):
                values[key] = os.environ[env]
        missing = [k for k in ("local_root", "cloud_remote") if not values.get(k)]
        if missing:
            raise ConfigError(
                f"cichlid-data is not configured on this machine (missing {', '.join(missing)}).\n"
                f"Run `cichlid-data init` to write {p} with the lab defaults, or set "
                f"CICHLID_DATA_ROOT and CICHLID_DATA_REMOTE.")
        return cls(local_root=Path(values["local_root"]).expanduser(), cloud_remote=values["cloud_remote"].rstrip("/"),
                   rclone=values.get("rclone", "rclone"))


def write_default_config(path=None, local_root=None, cloud_remote=None, overwrite=False):
    p = Path(path) if path else config_path()
    if p.exists() and not overwrite:
        raise ConfigError(f"{p} already exists (use --force to replace it)")
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(f'local_root = "{local_root or default_local_root()}"\n'
                 f'cloud_remote = "{cloud_remote or DEFAULT_REMOTE}"\n'
                 f'rclone = "rclone"\n')
    return p
