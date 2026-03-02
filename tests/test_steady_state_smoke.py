from __future__ import annotations

import sys
from pathlib import Path

import pyomo.environ as pyo

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from asa_cm_control.models.steady_state import build_model


def test_build_model_has_zero_dof() -> None:
    model = build_model()

    assert model.fs.feed_rate.fixed
    assert model.fs.conversion.fixed
    assert pyo.value(model.fs.product_rate) >= 0
