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


def resolve_extends(manifest: Mapping) -> dict:
    """Return the manifest with every case's "extends" parent merged in.

    A case with "extends": "<parent name>" is the parent (itself resolved first) overridden by the
    case's own keys: dicts merge recursively, every other value (lists included) replaces the
    parent's. The "extends" key is dropped. Cycles and missing parents raise ValueError; a
    manifest without "extends" comes back equal to the input.
    """
    cases = {case["name"]: case for case in manifest.get("cases", [])}

    def resolve(name: str, child: str | None, chain: tuple[str, ...]) -> dict:
        if name not in cases:
            raise ValueError(f"case {child!r} extends missing case {name!r}")
        if name in chain:
            raise ValueError("extends cycle: " + " -> ".join(chain + (name,)))
        case = cases[name]
        if "extends" not in case:
            return deepcopy(dict(case))
        own = {key: value for key, value in case.items() if key != "extends"}
        return _merge_dicts(resolve(case["extends"], name, chain + (name,)), own)

    resolved = dict(manifest)
    if "cases" in manifest:
        resolved["cases"] = [resolve(case["name"], None, ()) for case in manifest["cases"]]
    return resolved


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
