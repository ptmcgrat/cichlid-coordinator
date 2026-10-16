"""All transfers between this machine and the cloud.

    S = get_store()
    S.ensure(L.proteome("Mzeb"))          # fetch if missing locally; returns the local path
    S.pull("families/TreeFam9/*")         # unit names, or globs over them
    S.push("results/Compara2026/orthology")
    S.status()

Each unit's content fingerprint (paths, sizes, modification times) is recorded in the cloud index
(<cloud>/.cichlid-data/index.json) and in this machine's state (<local>/.cichlid-data/state.json):
    pull  downloads a unit if the cloud copy differs and the local copy is missing or unchanged
          since the last sync; if both changed it refuses (conflict) unless force
    push  uploads a unit if it changed locally and the cloud copy has not changed since the last
          sync; otherwise it refuses (conflict) unless force. Immutable units are never replaced
          by different content unless force. A push holds a cloud lock for the unit.
Every transfer is appended to <local>/.cichlid-data/transfers.log.
"""
import fnmatch
import hashlib
import json
import logging
import os
import shutil
import subprocess
import tarfile
import tempfile
import time
from pathlib import Path

from .provenance import stamp
from .units import Registry

logger = logging.getLogger(__name__)
INDEX = ".cichlid-data/index.json"
LOCKS = ".cichlid-data/locks"
LOCK_MAX_AGE = 6 * 3600
RCLONE_NOT_FOUND = frozenset({3, 4})
RCLONE_FATAL = frozenset({1, 7})


class TransferError(IOError):
    pass


# =============================================================================== backends
class LocalBackend:
    """A folder standing in for the cloud (tests, or a shared filesystem)."""

    def __init__(self, root):
        self.root = Path(root)

    def _p(self, rel):
        return self.root / rel.strip("/")

    def list(self, rel):
        d = self._p(rel)
        if not d.is_dir():
            return []
        return sorted(e.name + ("/" if e.is_dir() else "") for e in d.iterdir())

    def exists(self, rel):
        return self._p(rel).exists()

    def read_text(self, rel):
        p = self._p(rel)
        return p.read_text() if p.is_file() else None

    def write_text(self, rel, text):
        p = self._p(rel)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)

    def upload(self, local, rel_dir, immutable=False):
        """Copy a local file into a cloud folder, or a local folder's contents into it. With
        immutable, refuse (before copying anything) if an existing cloud file would change."""
        dest = self._p(rel_dir)
        local = Path(local)
        pairs = ([(f, dest / f.relative_to(local)) for f in local.rglob("*") if f.is_file()] if local.is_dir()
                 else [(local, dest / local.name)])
        if immutable:
            changed = [str(t.relative_to(self.root)) for f, t in pairs if t.exists() and t.read_bytes() != f.read_bytes()]
            if changed:
                raise TransferError(f"immutable: would modify {', '.join(changed[:5])}")
        dest.mkdir(parents=True, exist_ok=True)
        for f, t in pairs:
            t.parent.mkdir(parents=True, exist_ok=True)
            if not (immutable and t.exists()):
                shutil.copy2(f, t)

    def download(self, rel, local_dir):
        """Copy a cloud file into a local folder, or a cloud folder's contents into it."""
        src, dest = self._p(rel), Path(local_dir)
        dest.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            shutil.copytree(src, dest, dirs_exist_ok=True)
        elif src.is_file():
            shutil.copy2(src, dest / src.name)
        else:
            raise TransferError(f"not in the cloud: {rel}")

    def delete(self, rel):
        p = self._p(rel)
        if p.is_dir():
            shutil.rmtree(p)
        elif p.exists():
            p.unlink()

    def move(self, src, dst):
        s, d = self._p(src), self._p(dst)
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(s), str(d))


class RcloneBackend:
    """The lab Dropbox (or any rclone remote). Retries transient failures; not-found and fatal
    exit codes stop immediately."""

    def __init__(self, remote, rclone="rclone", retries=3, backoff=2.0):
        self.remote, self.rclone, self.retries, self.backoff = remote.rstrip("/"), rclone, retries, backoff

    def _u(self, rel):
        return f"{self.remote}/{rel.strip('/')}" if rel.strip("/") else self.remote

    def run(self, args, allow=(), stdin=None):
        proc = None
        for attempt in range(1, self.retries + 1):
            proc = subprocess.run([self.rclone] + list(args), capture_output=True, encoding="utf-8", input=stdin)
            if proc.returncode == 0 or proc.returncode in allow:
                return proc
            if proc.returncode in RCLONE_NOT_FOUND or proc.returncode in RCLONE_FATAL:
                break
            logger.warning("rclone %s failed (exit %s), attempt %d/%d: %s", args[0], proc.returncode, attempt,
                           self.retries, proc.stderr.strip())
            if attempt < self.retries:
                time.sleep(self.backoff * attempt)
        return proc

    def _ok(self, proc, what):
        if proc.returncode != 0 and ".DS_Store" not in (proc.stderr or ""):
            raise TransferError(f"{what} failed (rclone exit {proc.returncode}): {proc.stderr.strip()}")

    def list(self, rel):
        r = self.run(["lsf", self._u(rel)], allow=RCLONE_NOT_FOUND)
        if r.returncode in RCLONE_NOT_FOUND:
            return []
        self._ok(r, f"listing {rel}")
        return [x for x in r.stdout.split("\n") if x]

    def exists(self, rel):
        parent, _, name = rel.strip("/").rpartition("/")
        return any(x.rstrip("/") == name for x in self.list(parent))

    def read_text(self, rel):
        r = self.run(["cat", self._u(rel)], allow=RCLONE_NOT_FOUND)
        return r.stdout if r.returncode == 0 else None

    def write_text(self, rel, text):
        self._ok(self.run(["rcat", self._u(rel)], stdin=text), f"writing {rel}")

    def upload(self, local, rel_dir, immutable=False):
        flags = ["--multi-thread-streams", "16", "--multi-thread-cutoff", "100Mi"] + (["--immutable"] if immutable else [])
        r = self.run(["copy"] + flags + [str(local), self._u(rel_dir)])
        if immutable and r.returncode != 0 and "immutable" in (r.stderr or "").lower():
            raise TransferError(f"immutable: would modify files under {rel_dir}: {r.stderr.strip()[-300:]}")
        self._ok(r, f"upload of {local}")

    def download(self, rel, local_dir):
        """rclone copy puts a file into the folder, or a folder's contents into it (as LocalBackend)."""
        Path(local_dir).mkdir(parents=True, exist_ok=True)
        flags = ["--multi-thread-streams", "16", "--multi-thread-cutoff", "100Mi"]
        if not self.exists(rel):
            raise TransferError(f"not in the cloud: {rel}")
        self._ok(self.run(["copy"] + flags + [self._u(rel), str(local_dir)]), f"download of {rel}")

    def delete(self, rel):
        parent, _, name = rel.strip("/").rpartition("/")
        is_dir = (name + "/") in self.list(parent)
        self._ok(self.run(["purge" if is_dir else "deletefile", self._u(rel)]), f"deleting {rel}")

    def move(self, src, dst):
        """Server-side move. Dropbox refuses to move a folder holding too many files in one
        operation; then the folder's files are moved individually (still server-side, in
        parallel) and the emptied folders removed."""
        r = self.run(["moveto", self._u(src), self._u(dst)])
        if r.returncode == 0:
            return
        parent, _, name = src.strip("/").rpartition("/")
        if not self.exists(src) and self.exists(dst):            # Dropbox finished it despite the error
            logger.warning("move of %s reported an error but completed", src)
            return
        if (name + "/") not in self.list(parent):
            self._ok(r, f"moving {src}")                        # a file: nothing else to try
        logger.warning("whole-folder move of %s failed (%s); moving its files individually", src,
                       (r.stderr or "").strip().splitlines()[-1:] or "")
        self._ok(self.run(["move", self._u(src), self._u(dst), "--transfers", "32", "--checkers", "32",
                           "--delete-empty-src-dirs"]), f"moving the files of {src}")
        if (name + "/") in self.list(parent):                   # rclone may leave the top folder
            self.run(["rmdir", self._u(src)])


# =============================================================================== helpers
def _excluded(rel, excludes):
    return any(fnmatch.fnmatch(rel, p) or fnmatch.fnmatch(os.path.basename(rel), p) for p in excludes)


def unit_files(base, members, excludes=()):
    """Files (relative to base) under the members, sorted, minus excludes."""
    out = []
    for m in members:
        top = base / m
        if top.is_file():
            out.append(m)
            continue
        for dirpath, dirnames, files in os.walk(top):
            dirnames.sort()
            for f in sorted(files):
                rel = os.path.relpath(os.path.join(dirpath, f), base)
                if f != ".DS_Store" and not _excluded(rel, excludes):
                    out.append(rel)
    return out


def fingerprint(base, members, excludes=()):
    """Hash of (path, size, mtime) of every file; None if there are none."""
    files = unit_files(base, members, excludes)
    if not files:
        return None
    h = hashlib.sha256()
    for rel in files:
        st = (base / rel).stat()
        h.update(f"{rel}\0{st.st_size}\0{int(st.st_mtime)}\n".encode())
    return h.hexdigest()[:20]


def make_tar(tar_path, base, members, excludes=()):
    def keep(info):
        if os.path.basename(info.name) == ".DS_Store" or (info.isfile() and _excluded(info.name, excludes)):
            return None
        return info
    Path(tar_path).parent.mkdir(parents=True, exist_ok=True)
    tmp = str(tar_path) + ".part"
    with tarfile.open(tmp, "w") as tf:
        for m in members:
            tf.add(os.path.join(base, m), arcname=m, filter=keep)
    os.replace(tmp, tar_path)


def extract_tar(tar_path, dest):
    with tarfile.open(tar_path) as tf:
        real = os.path.realpath(dest)
        for m in tf.getmembers():
            target = os.path.realpath(os.path.join(dest, m.name))
            if target != real and not target.startswith(real + os.sep):
                raise TransferError(f"unsafe path in {tar_path}: {m.name}")
        try:
            tf.extractall(dest, filter="data")
        except TypeError:
            tf.extractall(dest)


# =============================================================================== store
class Store:
    def __init__(self, layout, backend, force=False, dry_run=False, context="", repo=None):
        self.layout, self.backend = layout, backend
        self.force, self.dry, self.context, self.repo = force, dry_run, context, repo
        self.state_path = layout.state_dir() / "state.json"
        self.state = json.loads(self.state_path.read_text()) if self.state_path.is_file() else {}
        self._index = None
        self._registry = None

    # ---- bookkeeping
    def index(self):
        if self._index is None:
            text = self.backend.read_text(INDEX)
            self._index = json.loads(text) if text and text.strip() else {}
        return self._index

    def registry(self):
        if self._registry is None:
            self._registry = Registry(self.layout, self.index())
        return self._registry

    def _save(self):
        if self.dry:
            return
        self.layout.state_dir().mkdir(parents=True, exist_ok=True)
        self.state_path.write_text(json.dumps(self.state, indent=2, sort_keys=True))
        if self._index is not None:
            self.backend.write_text(INDEX, json.dumps(self._index, indent=2, sort_keys=True))

    def _log(self, direction, key, result):
        if self.dry:
            return
        self.layout.state_dir().mkdir(parents=True, exist_ok=True)
        s = stamp(self.context)
        with open(self.layout.state_dir() / "transfers.log", "a") as fh:
            fh.write(f"{s['time']}\t{s['user']}@{s['machine']}\t{direction}\t{key}\t{result}\n")

    # ---- unit members
    def _units(self, names, bulk=False):
        if isinstance(names, str):
            names = [names]
        out = []
        for n in names:
            found = self.registry().matching(n)
            if not found:
                raise KeyError(f"no unit matches {n!r}; see `cichlid-data units`")
            out += [u for u in found if not (bulk and u.policy == "on_demand" and u.name != n)]
        return out

    def _shards(self, unit):
        if unit.packing != "shards":
            return [None]
        keys = set()
        d = self.layout.local_root / unit.rel
        if d.is_dir():
            keys |= {unit.shard_by(e.name, e.is_dir()) for e in d.iterdir() if not e.name.startswith(".")}
        keys |= {k[len(unit.name) + 1:] for k in self.index() if k.startswith(unit.name + "/")}
        return sorted(keys)

    def _members(self, unit, shard):
        """Paths relative to the unit's parent folder that make up the unit (or shard) locally."""
        base = self.layout.local_root / unit.parent
        if unit.packing == "shards":
            d = self.layout.local_root / unit.rel
            if not d.is_dir():
                return []
            return [f"{unit.base}/{e.name}" for e in sorted(d.iterdir())
                    if not e.name.startswith(".") and unit.shard_by(e.name, e.is_dir()) == shard]
        return [unit.base] if (base / unit.base).exists() else []

    def local_fp(self, unit, shard=None):
        m = self._members(unit, shard)
        return fingerprint(self.layout.local_root / unit.parent, m, unit.excludes) if m else None

    # ---- locks
    def _lock(self, key):
        path = f"{LOCKS}/{key.replace('/', '__')}.lock"
        text = self.backend.read_text(path)
        if text:
            info = json.loads(text)
            age = time.time() - info.get("epoch", 0)
            if age < LOCK_MAX_AGE and not self.force:
                raise TransferError(f"{key} is locked by {info.get('user')}@{info.get('machine')} since "
                                    f"{info.get('time')}; wait, or use --force if that push died")
        s = stamp(self.context)
        s["epoch"] = time.time()
        self.backend.write_text(path, json.dumps(s))
        return path

    # ---- pull / push one unit (or shard)
    def _pull_one(self, unit, shard):
        key = unit.key(shard)
        cloud = self.index().get(key)
        if not cloud:
            return "not in cloud"
        local, last = self.local_fp(unit, shard), self.state.get(key)
        if local == cloud["fingerprint"]:
            self.state[key] = local
            return "up to date"
        if local is not None and local != last and not self.force:
            return "CONFLICT (local changed since last sync; use --force to overwrite local)"
        if self.dry:
            return "would download"
        parent = self.layout.local_root / unit.parent
        parent.mkdir(parents=True, exist_ok=True)
        if unit.packing == "files":
            self.backend.download(unit.rel, parent / unit.base if self._cloud_is_dir(unit) else parent)
        else:
            with tempfile.TemporaryDirectory() as tmp:
                self.backend.download(unit.cloud_tar(shard), tmp)
                tar = Path(tmp) / Path(unit.cloud_tar(shard)).name
                extract_tar(tar, parent)
        self.state[key] = self.local_fp(unit, shard)
        return "downloaded"

    def _cloud_is_dir(self, unit):
        parent, _, name = unit.rel.rpartition("/")
        return (name + "/") in self.backend.list(parent)

    def _push_one(self, unit, shard):
        key = unit.key(shard)
        members = self._members(unit, shard)
        local = self.local_fp(unit, shard) if members else None
        if local is None:
            return "nothing local"
        cloud = (self.index().get(key) or {}).get("fingerprint")
        last = self.state.get(key)
        if local == cloud:
            self.state[key] = local
            return "up to date"
        additions_only = unit.policy == "immutable" and unit.packing == "files" and not self.force
        if cloud is not None and unit.policy == "immutable" and unit.packing != "files" and not self.force:
            return "REFUSED (immutable: the cloud already has a different version; use --force to replace it)"
        if cloud is not None and cloud != last and not self.force and not additions_only:
            return "CONFLICT (cloud changed since last sync; pull first, or --force to overwrite the cloud)"
        if self.dry:
            return "would upload"
        lock = self._lock(key)
        try:
            parent = self.layout.local_root / unit.parent
            if unit.packing == "files":
                src = parent / unit.base
                try:
                    self.backend.upload(src, unit.rel if src.is_dir() else unit.parent, immutable=additions_only)
                except TransferError as e:
                    if additions_only and "immutable" in str(e):
                        return f"REFUSED ({e}; existing files of an immutable unit can't change; use --force)"
                    raise
            else:
                with tempfile.TemporaryDirectory() as tmp:
                    tar = Path(tmp) / Path(unit.cloud_tar(shard)).name
                    make_tar(tar, parent, members, unit.excludes)
                    self.backend.upload(tar, str(Path(unit.cloud_tar(shard)).parent))
            self.index()[key] = {"fingerprint": local, **stamp(self.context, self.repo)}
            self.state[key] = local
        finally:
            self.backend.delete(lock)
        return "uploaded"

    # ---- public
    def pull(self, names, bulk=True):
        results = {}
        for u in self._units(names, bulk):
            for shard in self._shards(u):
                res = self._pull_one(u, shard)
                results[u.key(shard)] = res
                self._log("pull", u.key(shard), res)
        self._save()
        return results

    def push(self, names):
        results = {}
        for u in self._units(names):
            for shard in self._shards(u):
                res = self._push_one(u, shard)
                results[u.key(shard)] = res
                if res != "nothing local":
                    self._log("push", u.key(shard), res)
        self._save()
        bad = {k: v for k, v in results.items() if v.startswith(("CONFLICT", "REFUSED"))}
        if bad:
            raise TransferError("push incomplete; nothing was overwritten for: " +
                                "; ".join(f"{k}: {v}" for k, v in bad.items()))
        return results

    def ensure(self, path):
        """Make a path present locally (pulling the unit or shard that owns it). -> local Path."""
        rel = path.rel if hasattr(path, "rel") else str(Path(path).relative_to(self.layout.local_root))
        local = self.layout.local_root / rel
        if local.exists() and not self.force:
            return local
        unit, shard = self.registry().for_path(rel)
        if unit.packing == "files" and rel != unit.rel:
            if not self.dry:                       # fetch just this file or folder
                self.backend.download(rel, local.parent if not self._is_cloud_dir(rel) else local)
                self._log("pull", rel, "downloaded (single path)")
            return local
        res = self._pull_one(unit, shard)
        self._log("pull", unit.key(shard), res)
        self._save()
        if not local.exists() and not self.dry:
            raise TransferError(f"{rel} is not in the cloud copy of {unit.key(shard)} ({res})")
        return local

    def _is_cloud_dir(self, rel):
        parent, _, name = rel.strip("/").rpartition("/")
        return (name + "/") in self.backend.list(parent)

    def status(self, names="*"):
        """-> list of (key, status, last uploader)"""
        rows = []
        for u in self._units(names):
            for shard in self._shards(u):
                key = u.key(shard)
                local = self.local_fp(u, shard)
                info = self.index().get(key) or {}
                cloud = info.get("fingerprint")
                st = ("in sync" if local and local == cloud else "local only" if local and not cloud
                      else "cloud only" if cloud and not local
                      else ("local changed" if self.state.get(key) == cloud else
                            "cloud changed" if self.state.get(key) == local else "both changed") if local else "-")
                rows.append((key, st, f"{info.get('user', '')}@{info.get('machine', '')} {info.get('time', '')}" if info else ""))
        return rows


def get_store(settings=None, **kw):
    from .config import Settings
    from .layout import Layout
    s = settings or Settings.load()
    return Store(Layout(s.local_root, s.cloud_remote), RcloneBackend(s.cloud_remote, s.rclone), **kw)
