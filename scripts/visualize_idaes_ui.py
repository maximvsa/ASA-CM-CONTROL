from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from idaes.core.base.flowsheet_model import UI
from asa_cm_control.steady_state_demo import build_model


def main() -> int:
    parser = argparse.ArgumentParser(description="Launch IDAES built-in flowsheet visualizer.")
    parser.add_argument("--title", default="ASA CM Demo", help="Title shown in the visualizer.")
    parser.add_argument(
        "--save",
        default="asa_cm_demo_vis.json",
        help="Filename for saved visualization JSON.",
    )
    parser.add_argument(
        "--save-dir",
        default=str(REPO_ROOT / "artifacts" / "flowsheet"),
        help="Directory for saved visualization JSON.",
    )
    parser.add_argument(
        "--no-block",
        action="store_true",
        help="Exit immediately after launching (server will stop when process exits).",
    )
    args = parser.parse_args()

    if not UI().installed:
        print("idaes_ui is not installed in this environment.")
        print("Install it with: pip install idaes-ui")
        return 1

    model = build_model()
    result = model.fs.visualize(args.title, save=args.save, save_dir=args.save_dir)

    print("IDAES visualizer launched.")
    print(f"Saved JSON: {Path(args.save_dir) / args.save}")
    if getattr(result, "port", None) is not None:
        print(f"Open in browser: http://localhost:{result.port}")
    print(f"Visualizer result: {result}")

    if args.no_block:
        return 0

    print("Visualizer server is running. Press Ctrl+C to stop.")
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("Stopping visualizer server.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
