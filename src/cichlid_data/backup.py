"""Helpers for scripts that write into the data tree: push what they wrote, and say where it went.

    from cichlid_data.backup import backup
    print(f"Backed up to {backup(out)}")      # pushes the unit owning `out`; returns its cloud path
"""
from pathlib import Path

from .config import Settings
from .layout import get_layout
from .store import get_store


def settings():
    return Settings.load()


def _rel(layout, path):
    return str(Path(path).resolve().relative_to(Path(layout.local_root).resolve()))


def cloud_path(path, layout=None):
    layout = layout or get_layout()
    return layout.path(_rel(layout, path)).cloud


def backup(path, store=None):
    """Push the unit (or shard) owning `path`. -> the path's cloud location."""
    store = store or get_store(context="backup")
    rel = _rel(store.layout, path)
    unit, _shard = store.registry().for_path(rel)
    store.push(unit.name)
    return store.layout.path(rel).cloud


