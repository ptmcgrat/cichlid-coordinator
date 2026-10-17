"""The documents say what the code does: the specification's command list against the real command line."""
import io
import re
import unittest
from contextlib import redirect_stdout

from cichlid_data import __version__
from cichlid_data.cli import PASSTHROUGH, main as cli
from cichlid_data.config import config_dir

SPEC = config_dir().parent / "docs" / "coordinator_spec.md"


def spec_lines():
    block = re.search(r"## 5\. Command line\n\n```\n(.*?)```", SPEC.read_text(), re.S)
    return [ln for ln in block.group(1).splitlines() if ln.strip()]


def spec_verbs(lines):
    """cichlid-data units | status ...  has two verbs; cichlid-data provenance write | check ... has one."""
    out = {}
    for ln in lines:
        t = ln.split()
        verbs = [t[1]] + ([t[3]] if len(t) > 3 and t[2] == "|" else [])
        for v in verbs:
            out.setdefault(v, []).append(ln)
    return out


def help_text(*argv):
    f = io.StringIO()
    with redirect_stdout(f):
        try:
            cli(list(argv))
        except SystemExit:
            pass
    return f.getvalue()


@unittest.skipUnless(SPEC.is_file(), "docs/coordinator_spec.md not present")
class TestSpecCommands(unittest.TestCase):
    def test_every_command_is_documented_and_nothing_else(self):
        usage = re.search(r"\{([a-z\-,]+)\}", help_text("--help")).group(1).split(",")
        code = set(usage) | {v for v, _ in PASSTHROUGH}
        self.assertGreaterEqual(len(code), 20)
        self.assertEqual(code, set(spec_verbs(spec_lines())))

    def test_every_option_is_documented(self):
        verbs = spec_verbs(spec_lines())
        for verb in re.search(r"\{([a-z\-,]+)\}", help_text("--help")).group(1).split(","):
            if verb == "prepare-local":                         # its options are the preparation tool's own
                continue
            options = set(re.findall(r"(--[a-z][a-z0-9\-]*)", help_text(verb, "--help"))) - {"--help"}
            text = " ".join(verbs[verb])
            self.assertEqual([o for o in sorted(options) if o not in text], [], f"`{verb}` options missing from the spec")



class TestVersion(unittest.TestCase):
    def test_pyproject_and_package_agree(self):
        """HANDOFF.md: bump pyproject.toml and src/cichlid_data/__init__.py together."""
        pyproject = SPEC.parents[1] / "pyproject.toml"
        if not pyproject.is_file():
            self.skipTest("pyproject.toml not present")
        self.assertIn(f'version = "{__version__}"', pyproject.read_text())


if __name__ == "__main__":
    unittest.main()
