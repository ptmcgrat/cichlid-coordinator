"""File locations, cloud transfers and inputs for the cichlid protein phylogeny project.

Public API: Settings, Layout, DataPath, ResultsLayout, DataUnit, Registry, Store, Catalog, Species,
get_layout, get_store.
"""
__version__ = "0.12.0"

from .config import ConfigError, Settings  # noqa: E402
from .layout import DataPath, Layout, ResultsLayout, get_layout  # noqa: E402
from .clades import Clades  # noqa: E402
from . import master_tree  # noqa: E402,F401
from .store import LocalBackend, RcloneBackend, Store, TransferError, get_store  # noqa: E402
from .units import DataUnit, Registry  # noqa: E402
from .catalog import Catalog, Species  # noqa: E402
from . import schemas  # noqa: E402

__all__ = ["Settings", "ConfigError", "Layout", "DataPath", "ResultsLayout", "get_layout", "DataUnit", "Registry",
           "Store", "LocalBackend", "RcloneBackend", "TransferError", "get_store", "Catalog", "Species", "Clades", "schemas",
           "__version__"]


