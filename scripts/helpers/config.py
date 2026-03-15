import sys
from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass
class Style:
    dpi: int
    axis_label_fontsize: int
    axis_label_fontweight: str
    tick_fontsize: int
    tick_fontweight: str
    kpam_color: str
    lit_color: str
    pyruvylation_color: str
    acetylation_color: str
    gray_color: str
    both_color: str


class Config:
    _CONFIG_PATH = Path(__file__).resolve().parents[2] / "config" / "config.yml"
    _ANALYSIS_SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "analysis"

    def __init__(self):
        with open(self._CONFIG_PATH) as f:
            data = yaml.safe_load(f)
        paths = data["paths"]
        self.input_dir = Path(paths["input_dir"])
        self.output_dir = Path(paths["output_dir"])

        s = data.get("style", {})
        self.style = Style(
            dpi=s.get("dpi", 200),
            axis_label_fontsize=s.get("axis_label_fontsize", 12),
            axis_label_fontweight=s.get("axis_label_fontweight", "bold"),
            tick_fontsize=s.get("tick_fontsize", 8),
            tick_fontweight=s.get("tick_fontweight", "bold"),
            kpam_color=s.get("kpam_color", "#1f77b4"),
            lit_color=s.get("lit_color", "#d62728"),
            pyruvylation_color=s.get("pyruvylation_color", "#4f81bd"),
            acetylation_color=s.get("acetylation_color", "#f2c94c"),
            gray_color=s.get("gray_color", "#bfbfbf"),
            both_color=s.get("both_color", "#DC143C"),
        )

        a = data.get("analysis", {})
        self.network_cutoff = float(a.get("network_cutoff", 0.5))

        self._check_analysis_readmes()

    def _check_analysis_readmes(self):
        if not self._ANALYSIS_SCRIPTS_DIR.exists():
            return
        for chapter_dir in sorted(self._ANALYSIS_SCRIPTS_DIR.iterdir()):
            if not chapter_dir.is_dir() or chapter_dir.name.startswith("."):
                continue
            if not (chapter_dir / "main.py").exists():
                continue
            readme = chapter_dir / "README.md"
            if not readme.exists():
                print(
                    f"ERROR: Missing README.md in scripts/analysis/{chapter_dir.name}/\n"
                    f"Every analysis chapter folder must contain a README.md describing "
                    f"the input data, output tables, and figures.",
                    file=sys.stderr,
                )
                sys.exit(1)
