#!/usr/bin/env python3
"""A stand-in for rclone over a local folder: remote paths 'fake:/x' map to $FAKE_CLOUD/x."""
import os
import shutil
import sys
from pathlib import Path

root = Path(os.environ["FAKE_CLOUD"])


def p(arg):
    return root / arg.split(":", 1)[1].lstrip("/") if ":" in arg and not arg.startswith("/") else Path(arg)


args = [a for a in sys.argv[1:]]
cmd = args[0]
rest, i = [], 1
while i < len(args):                                   # drop flags and their values
    if args[i].startswith("--"):
        i += 2 if args[i] in ("--multi-thread-streams", "--multi-thread-cutoff") else 1
        continue
    rest.append(args[i])
    i += 1
if cmd == "lsf":
    d = p(rest[0])
    if not d.is_dir():
        sys.exit(3)
    for e in sorted(d.iterdir()):
        print(e.name + ("/" if e.is_dir() else ""))
elif cmd == "cat":
    f = p(rest[0])
    if not f.is_file():
        sys.exit(3)
    sys.stdout.write(f.read_text())
elif cmd == "rcat":
    f = p(rest[0])
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(sys.stdin.read())
elif cmd == "copy":
    src, dst = p(rest[0]), p(rest[1])
    if "--immutable" in sys.argv:                              # as rclone: existing files must not change
        pairs = [(x, dst / x.relative_to(src)) for x in src.rglob("*") if x.is_file()] if src.is_dir() else [(src, dst / src.name)]
        bad = [t for s_, t in pairs if t.exists() and t.read_bytes() != s_.read_bytes()]
        if bad:
            sys.stderr.write(f"ERROR : {bad[0].name}: Source and destination exist but do not match: immutable file modified\n")
            sys.exit(1)
    dst.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True)
    elif src.is_file():
        shutil.copy2(src, dst / src.name)
    else:
        sys.exit(3)
elif cmd == "moveto":
    src = p(rest[0])
    if src.is_dir() and os.environ.get("FAKE_FAIL_MOVEDIR"):     # like Dropbox with too many files
        sys.stderr.write("ERROR : Server side directory move failed: MoveDir failed:\n")
        sys.exit(1)
    d = p(rest[1])
    d.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(d))
    if os.environ.get("FAKE_MOVE_ERRS_BUT_COMPLETES"):          # as Dropbox did with NewTreeFiles
        sys.stderr.write("ERROR : Server side directory move failed: MoveDir failed:\n")
        sys.exit(1)
elif cmd == "move":                                            # file by file, then empty folders
    src, dst = p(rest[0]), p(rest[1])
    for f in sorted(x for x in src.rglob("*") if x.is_file()):
        t = dst / f.relative_to(src)
        t.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(f), str(t))
    if "--delete-empty-src-dirs" in sys.argv:
        for d in sorted((x for x in src.rglob("*") if x.is_dir()), key=lambda x: -len(x.parts)):
            d.rmdir()
elif cmd == "rmdir":
    p(rest[0]).rmdir()
elif cmd == "deletefile":
    p(rest[0]).unlink()
elif cmd == "purge":
    shutil.rmtree(p(rest[0]))
else:
    sys.exit(1)


