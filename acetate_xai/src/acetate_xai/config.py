from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when a config file cannot be loaded or has invalid structure."""


def load_config(path: str | Path) -> dict[str, Any]:
    """
    Load a YAML or JSON config file into a dict.

    Parameters
    ----------
    path:
        Path to a .yaml/.yml or .json file.
    """
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"Config file not found: {p}")

    suffix = p.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        with p.open("r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
    elif suffix == ".json":
        with p.open("r", encoding="utf-8") as f:
            data = json.load(f)
    else:
        raise ConfigError(f"Unsupported config extension: {suffix} (expected .yaml/.yml/.json)")

    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"Config must be a mapping/dict, got: {type(data).__name__}")
    return data


@dataclass(frozen=True)
class Paths:
    """Convenience container for project paths."""

    project_root: Path

    @property
    def configs_dir(self) -> Path:
        return self.project_root / "configs"

    @property
    def data_dir(self) -> Path:
        return self.project_root / "data"

    @property
    def results_dir(self) -> Path:
        return self.project_root / "results"

