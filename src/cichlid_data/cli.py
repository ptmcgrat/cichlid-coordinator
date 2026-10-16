"""cichlid-data: file locations, transfers and the layout migration.

    cichlid-data init [--root DIR] [--remote REMOTE] [--force]
    cichlid-data where proteome Mzeb          # print a local path (--cloud for the cloud path)
    cichlid-data env                          # shell exports:  eval "$(cichlid-data env)"
    cichlid-data units [PATTERN]
    cichlid-data status [PATTERN]
    cichlid-data pull PATTERN... [--force] [--dry-run]
    cichlid-data push PATTERN... [--force] [--dry-run]
    cichlid-data ensure RELATIVE_PATH
    cichlid-data retire UNIT_KEY                     # delete a unit's or shard's cloud copy
    cichlid-data download genomes [--species P,..] [--check-only] [--push]   # Stage 1: NCBI or local files
    cichlid-data prepare-local --prefix P --cds ... # convert a non-NCBI annotation (Etroplus)
    cichlid-data download annotations|expression|populations [options]
    cichlid-data populations samples [options]       # sample metadata for the population VCFs
    cichlid-data prepare compara [options]           # prune the Compara dumps to the reference panel
    cichlid-data summarize h5ad|rnaseq [options]     # expression summaries
      (these run the moved scripts with their own options: add --help to see them)
    cichlid-data proteomes [--species P,..] [--push] # one protein per gene, for each species
    cichlid-data check [--method M] [--fetch] [--deep]   # is everything each method needs here and sound?
    cichlid-data validate PATH... [--kind K]        # check files against their schemas
    cichlid-data catalog build|show|missing|validate # the species catalog (--method, --role)
    cichlid-data provenance write|check|show        # source records for every input
    cichlid-data migrate [--apply] [--cleanup]       # moves the inputs to the new layout

PATTERN is a unit name or a glob over unit names, e.g. "results/Compara2026/*".
"""
import argparse
import sys

from . import __version__
from .config import ConfigError, Settings, write_default_config
from .layout import Layout
from .migrate import Migration
from .store import RcloneBackend, Store, TransferError

WHERE = {
    "root": lambda L, a: L.path(""),
    "proteome": lambda L, a: L.proteome(a[0]),
    "genes": lambda L, a: L.genes_table(a[0]),
    "download": lambda L, a: L.download(a[0]),
    "local-source": lambda L, a: L.local_source(a[0]),
    "proteomes": lambda L, a: L.proteomes(),
    "catalog": lambda L, a: L.catalog(),
    "family-resource": lambda L, a: L.family_resource(a[0]),
    "annotations": lambda L, a: L.annotations(a[0] if a else ""),
    "expression": lambda L, a: L.expression(a[0] if a else ""),
    "population": lambda L, a: L.population(a[0] if a else ""),
    "results": lambda L, a: L.results(a[0]).part(a[1]) if len(a) > 1 else L.results(a[0]).root(),
    "analysis": lambda L, a: L.analysis(a[0]),
    "archive": lambda L, a: L.archive(a[0] if a else ""),
}


def build(args, settings):
    L = Layout(settings.local_root, settings.cloud_remote)
    backend = getattr(args, "backend", None) or RcloneBackend(settings.cloud_remote, settings.rclone)
    return L, backend


def main(argv=None, backend=None, settings=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv[:1] == ["prepare-local"]:          # every option belongs to the preparation tool
        return prepare_local(argv[1:], settings)
    if tuple(argv[:2]) in PASSTHROUGH:          # moved scripts keep their own options
        return passthrough(tuple(argv[:2]), argv[2:], settings)
    ap = argparse.ArgumentParser(prog="cichlid-data", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("init")
    p.add_argument("--root")
    p.add_argument("--remote")
    p.add_argument("--force", action="store_true")
    p = sub.add_parser("where")
    p.add_argument("kind", choices=sorted(WHERE))
    p.add_argument("args", nargs="*")
    p.add_argument("--cloud", action="store_true")
    sub.add_parser("env")
    p = sub.add_parser("units")
    p.add_argument("pattern", nargs="?", default="*")
    p = sub.add_parser("status")
    p.add_argument("pattern", nargs="?", default="*")
    for name in ("pull", "push"):
        p = sub.add_parser(name)
        p.add_argument("patterns", nargs="+")
        p.add_argument("--force", action="store_true")
        p.add_argument("--dry-run", action="store_true")
    p = sub.add_parser("retire", help="delete a unit's or shard's cloud copy and index entry")
    p.add_argument("key", help="e.g. proteins/downloads/_top")
    p.add_argument("--dry-run", action="store_true")
    p = sub.add_parser("ensure")
    p.add_argument("path")
    p = sub.add_parser("download", help="download inputs")
    p.add_argument("source", choices=["genomes"])
    p.add_argument("--species", help="comma-separated prefixes (default: every species in species.tsv)")
    p.add_argument("--check-only", action="store_true", help="check each species against NCBI; download nothing")
    p.add_argument("--force", action="store_true", help="download again even if present and verified")
    p.add_argument("--keep-zip", action="store_true")
    p.add_argument("--push", action="store_true", help="push the downloads to the cloud afterwards")
    p.add_argument("--no-fetch", action="store_true", help="don't pull species already in the cloud")
    sub.add_parser("prepare-local", help="convert a non-NCBI annotation (e.g. Etroplus); options as prepare_local_species")
    p = sub.add_parser("proteomes", help="one protein per gene, for each species")
    p.add_argument("--species", help="comma-separated prefixes (default: every species in species.tsv)")
    p.add_argument("--select", choices=["complete-longest", "longest"], default="complete-longest")
    p.add_argument("--force", action="store_true", help="rebuild even if up to date")
    p.add_argument("--push", action="store_true", help="push the proteomes to the cloud afterwards")
    p.add_argument("--no-fetch", action="store_true", help="don't pull missing downloads from the cloud")
    p = sub.add_parser("check", help="is everything each method needs here and sound?")
    p.add_argument("--method", help="one method (default: all)")
    p.add_argument("--fetch", action="store_true", help="pull anything missing from the cloud first")
    p.add_argument("--deep", action="store_true", help="also compare checksums with source records")
    p = sub.add_parser("validate", help="check files against their schemas")
    p.add_argument("paths", nargs="+")
    p.add_argument("--kind", choices=["genes", "assign", "orthologs", "proteome", "run", "family"])
    p.add_argument("--prefix")
    p = sub.add_parser("catalog", help="the species catalog")
    p.add_argument("action", choices=["build", "show", "missing", "validate"])
    p.add_argument("--method")
    p.add_argument("--role")
    p = sub.add_parser("provenance", help="source records for every input")
    p.add_argument("action", choices=["write", "check", "show"])
    p.add_argument("path", nargs="?", help="show: a folder relative to the data root")
    p.add_argument("--overwrite", action="store_true", help="write: replace existing records")
    p.add_argument("--md5-limit-gb", type=float, default=1.0, help="write: checksum files up to this size")
    p.add_argument("--deep", action="store_true", help="check: also compare checksums")
    p = sub.add_parser("migrate")
    p.add_argument("--apply", action="store_true")
    p.add_argument("--cleanup", action="store_true")
    args = ap.parse_args(argv)

    if args.cmd == "init":
        try:
            path = write_default_config(local_root=args.root, cloud_remote=args.remote, overwrite=args.force)
        except ConfigError as e:
            sys.exit(str(e))
        print(f"Wrote {path}")
        return 0
    try:
        settings = settings or Settings.load()
    except ConfigError as e:
        sys.exit(str(e))
    args.backend = backend
    L, B = build(args, settings)

    if args.cmd == "where":
        try:
            p = WHERE[args.kind](L, args.args)
        except IndexError:
            sys.exit(f"`where {args.kind}` needs an argument (e.g. a species prefix or name)")
        print(p.cloud if args.cloud else p.local)
        return 0
    if args.cmd == "download":
        from .catalog import Catalog
        from .config import config_dir
        from .downloaders import ncbi
        store = None if args.no_fetch else Store(L, B, context="download")
        try:
            summary, failed = ncbi.download(L, config_dir() / "species.tsv",
                                            species=args.species.split(",") if args.species else None,
                                            check_only=args.check_only, force=args.force, keep_zip=args.keep_zip,
                                            store=store)
        except ncbi.ManifestError as e:
            sys.exit(str(e))
        if not args.check_only:
            C = Catalog.build(L)
            if not C.validate()[0]:
                C.write(L)
                print(f"Catalog updated: {sum(s.status in ('downloaded', 'proteome ready') for s in C.organisms.values())} "
                      f"species downloaded")
            if args.push and any(r["status"] == "downloaded" for r in summary):
                res = Store(L, B, context="download").push("proteins/downloads")
                print("pushed proteins/downloads: " + ", ".join(sorted(set(res.values()))))
        if failed:
            sys.exit(f"{len(failed)} species failed: {', '.join(failed)}")
        return 0
    if args.cmd == "proteomes":
        from .catalog import Catalog
        from .config import config_dir
        from .proteomes import build as build_proteomes
        store = None if args.no_fetch else Store(L, B, context="proteomes")
        summaries, failed = build_proteomes(L, config_dir() / "species.tsv",
                                            species=args.species.split(",") if args.species else None,
                                            mode=args.select, force=args.force, store=store)
        C = Catalog.build(L)
        if not C.validate()[0]:
            C.write(L)
            print(f"Catalog updated: {sum(s.status == 'proteome ready' for s in C.organisms.values())} proteomes ready")
        if args.push and summaries:
            res = Store(L, B, context="proteomes").push("proteins/proteomes")
            print("pushed proteins/proteomes: " + ", ".join(sorted(set(res.values()))))
        if failed:
            sys.exit(f"{len(failed)} species failed: {', '.join(failed)}")
        return 0
    if args.cmd == "check":
        from .catalog import Catalog
        from .check import check as run_check
        C = Catalog.build(L)
        store = Store(L, B, context="check") if args.fetch else None
        findings = run_check(L, C, [args.method] if args.method else None, fetch=args.fetch, deep=args.deep,
                             store=store, say=print)
        bad = [f for f in findings if not f.ok]
        for m in dict.fromkeys(f.method for f in findings):
            mine = [f for f in findings if f.method == m]
            n_bad = sum(not f.ok for f in mine)
            print(f"{m}: {len(mine) - n_bad} of {len(mine)} inputs ready" + ("" if not n_bad else f"; {n_bad} problem(s):"))
            for f in mine:
                if not f.ok:
                    print(f"    {f.item}: {f.detail}")
        if bad:
            sys.exit(1)
        return 0
    if args.cmd == "validate":
        from .schemas import SchemaError, validate as _validate
        failed = 0
        for path in args.paths:
            try:
                problems = _validate(path, kind=args.kind, prefix=args.prefix)
            except (SchemaError, OSError) as e:
                problems = [str(e)]
            failed += bool(problems)
            print(f"{path}: " + ("ok" if not problems else f"{len(problems)} problem(s)"))
            for pr in problems[:10]:
                print(f"    {pr}")
        if failed:
            sys.exit(f"{failed} of {len(args.paths)} failed validation")
        return 0
    if args.cmd == "catalog":
        from .catalog import Catalog
        if args.action == "build":
            C = Catalog.build(L)
            errors, warnings = C.validate()
            for w in warnings:
                print(f"  warning: {w}")
            if errors:
                sys.exit("catalog not written:\n  " + "\n  ".join(errors))
            path = C.write(L)
            ready = sum(s.status == "proteome ready" for s in C.organisms.values())
            print(f"Wrote {path}: {len(C.organisms)} organisms ({ready} with proteomes ready); methods: " +
                  ", ".join(f"{n} {len(m['species'])}" for n, m in C.methods.items()))
            return 0
        C = Catalog.build(L) if args.action == "validate" else Catalog.load(L)
        if args.action == "validate":
            errors, warnings = C.validate()
            for w in warnings:
                print(f"  warning: {w}")
            for e in errors:
                print(f"  ERROR: {e}")
            print("catalog configuration is valid" if not errors else "")
            if errors:
                sys.exit(1)
        elif args.action == "missing":
            if not args.method:
                sys.exit("catalog missing needs --method, e.g. --method OrthoFinder2026")
            miss = C.missing(args.method)
            print(f"{args.method}: {len(C.species(method=args.method)) - len(miss)} of "
                  f"{len(C.species(method=args.method))} proteomes ready; missing:")
            for s in miss:
                print(f"  {s.prefix:<6} {s.species:<30} {s.status:<18} {'in species.tsv' if s.in_manifest else 'not in species.tsv'}")
            if not miss:
                print("  none")
        else:
            for s in C.species(method=args.method, role=args.role):
                print(f"  {s.prefix:<6} {s.species:<30} {s.role:<10} {s.status:<18} {s.accession:<18} "
                      f"{s.counts.get('proteins', '')!s:>7}  {','.join(s.used_by)}")
        return 0
    if args.cmd == "provenance":
        import json as _json
        from .provenance import read_record
        from .provenance_cmd import check_records, write_records
        if args.action == "write":
            write_records(L, overwrite=args.overwrite, md5_limit=int(args.md5_limit_gb * (1 << 30)))
        elif args.action == "check":
            if check_records(L, deep=args.deep):
                sys.exit("provenance check found problems")
        else:
            if not args.path:
                sys.exit("provenance show needs a folder, e.g. Inputs/Families/TreeFam9")
            rec = read_record(L.path(args.path).local)
            if rec is None:
                sys.exit(f"no source record in {args.path}")
            d = rec.to_dict()
            d["files"] = f"{len(rec.files)} files"
            print(_json.dumps(d, indent=2))
        return 0
    if args.cmd == "env":
        for k, v in L.env().items():
            print(f"export {k}='{v}'")
        return 0
    S = Store(L, B, force=getattr(args, "force", False), dry_run=getattr(args, "dry_run", False), context="cli")
    try:
        if args.cmd == "units":
            for u in S.registry().matching(args.pattern):
                print(f"{u.name:<42} {u.policy:<12} {u.packing:<7} {u.rel}")
        elif args.cmd == "status":
            print(f"local: {L.local_root}\ncloud: {L.cloud_remote}\n")
            for key, st, who in S.status(args.pattern):
                print(f"  {key:<46} {st:<14} {who}")
        elif args.cmd in ("pull", "push"):
            res = getattr(S, args.cmd)(args.patterns)
            for k, v in res.items():
                if v != "nothing local":
                    print(f"  {k}: {v}")
        elif args.cmd == "retire":
            print(f"{args.key}: {S.retire(args.key)}")
        elif args.cmd == "ensure":
            print(S.ensure(L.path(args.path)))
        elif args.cmd == "migrate":
            Migration(L, B).run(lambda: Store(L, B, context="migrate"), apply=args.apply, cleanup=args.cleanup)
    except (TransferError, KeyError) as e:
        sys.exit(str(e))
    return 0


# command -> (module, help). Their options are the moved scripts' own; run with --help to see them.
PASSTHROUGH = {
    ("download", "annotations"): ("annotations", "GO and PANTHER annotation files"),
    ("download", "expression"): ("expression", "HPA, Bgee or Daniocell expression data (--source)"),
    ("download", "populations"): ("blumer", "Blumer et al. 2025 population VCFs and tables"),
    ("populations", "samples"): ("blumer_samples", "sample metadata for the population VCFs"),
    ("prepare", "compara"): ("compara", "prune the Compara dumps to the reference panel"),
    ("summarize", "h5ad"): ("h5ad", "single-cell atlas (h5ad) to per-cell-type summaries; needs the [expression] extra"),
    ("summarize", "rnaseq"): ("rnaseq", "tissue RNA-seq (nf-core outputs) to per-tissue summaries"),
}


def passthrough(command, argv, settings=None):
    """Run a moved script with its own options, adding the coordinator's default paths it requires."""
    import importlib
    module = importlib.import_module(f".downloaders.{PASSTHROUGH[command][0]}", __package__)
    argv = list(argv)
    if command == ("prepare", "compara") and "--help" not in argv and "-h" not in argv:
        from .config import config_dir
        try:
            s = settings or Settings.load()
        except ConfigError as e:
            sys.exit(str(e))
        C = Layout(s.local_root, s.cloud_remote).family_resource("Compara116").local
        defaults = {"--nhx": C / "dumps" / "Compara.116.protein_default.nhx.emf.gz",
                    "--aln": C / "dumps" / "Compara.116.protein_default.aa.fasta.gz",
                    "--panel": config_dir() / "compara_panel.tsv", "--outdir": C / "prepared"}
        for opt, path in defaults.items():
            if opt not in argv:
                argv += [opt, str(path)]
    module.main(argv)
    return 0


def prepare_local(argv, settings=None):
    """Run the local-species preparation; --outdir defaults to the species' local-source folder."""
    from .downloaders import local
    argv = list(argv)
    if "--outdir" not in argv and "--prefix" in argv and argv.index("--prefix") + 1 < len(argv):
        try:
            s = settings or Settings.load()
        except ConfigError as e:
            sys.exit(str(e))
        prefix = argv[argv.index("--prefix") + 1]
        argv += ["--outdir", str(Layout(s.local_root, s.cloud_remote).local_source(prefix).local / "prepared")]
    local.main(argv)
    return 0


if __name__ == "__main__":
    main()
