__version__ = "0.1.0"

from .cli import main as hfr_main
from .compute_single import main as compute_folder_main
from .multi.hfr_driver import main_cli as multi_main
from .multi.run_multiple import main_cli as multi_run_main
from .multi.compute_multiple import main_cli as multi_compute_main

__all__ = [
    "__version__",
    "hfr_main",
    "compute_folder_main",
    "multi_main",
    "multi_run_main",
    "multi_compute_main",
]
