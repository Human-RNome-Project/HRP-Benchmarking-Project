"""Shared path definitions for the draft reference pipeline."""
from pathlib import Path

ROOT = Path(__file__).parent.parent
INPUTS = ROOT / "inputs"
OUTPUTS = ROOT / "outputs"
THRESHOLDS = OUTPUTS / "thresholds"
FILTERED = OUTPUTS / "filtered_platform_beds"
STATE = OUTPUTS / "state_assignments"
DRAFT = OUTPUTS / "draft_reference_beds"
TIERED = OUTPUTS / "tiered_lists"
TIERED_TRNA = OUTPUTS / "tiered_tRNA"
