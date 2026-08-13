"""Small manifest-level namelist default handling shared by test runners."""

from __future__ import annotations

from copy import deepcopy
from collections.abc import Mapping


def _merge_dicts(defaults: Mapping, overrides: Mapping) -> dict:
    """Recursively merge mapping values, with overrides taking precedence."""
    merged = deepcopy(dict(defaults))
    for key, value in overrides.items():
        if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
            merged[key] = _merge_dicts(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def apply_manifest_defaults(manifest: Mapping, case: Mapping) -> dict:
    """Return a case whose namelists include optional manifest-level defaults."""
    defaults = manifest.get("defaults", {})
    default_namelists = defaults.get("namelists", {}) if isinstance(defaults, Mapping) else {}
    case_namelists = case.get("namelists", {})
    if not isinstance(default_namelists, Mapping) or not isinstance(case_namelists, Mapping):
        raise TypeError("manifest defaults and case namelists must be mappings")

    merged_case = deepcopy(dict(case))
    merged_case["namelists"] = _merge_dicts(default_namelists, case_namelists)
    return merged_case
