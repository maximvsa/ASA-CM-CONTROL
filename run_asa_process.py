"""Root launcher for running the ASA process flowsheet.

Run from repository root:
    python run_asa_process.py
"""

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parent
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from asa_cm_control.asa_process_flowsheet import main


if __name__ == "__main__":
    main()
