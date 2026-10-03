"""Helpers for optional mature analysis backends used by immune."""

from __future__ import annotations

from importlib import import_module
from importlib.util import find_spec
from types import ModuleType


def require_dependency(module: str, *, extra: str, feature: str) -> ModuleType:
    """Import an optional dependency or raise an actionable error."""

    try:
        return import_module(module)
    except ImportError as error:
        raise ImportError(
            f"{feature} requires the optional dependency {module!r}. "
            f"Install it with `python -m pip install -e '.[{extra}]'`."
        ) from error


def backend_status() -> dict[str, bool]:
    """Report which supported Python analysis backends are available."""

    modules = {
        "scikit-bio": "skbio",
        "scirpy": "scirpy",
        "pydeseq2": "pydeseq2",
    }
    return {name: find_spec(module) is not None for name, module in modules.items()}
