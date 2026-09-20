"""Normalize the small public interface into existing stage-specific options."""
from dataclasses import dataclass
from collections.abc import Mapping
from os import PathLike
from pathlib import Path

from .discovery import kind_name
from .options import INPUT, OUTPUT, RUNTIME_DEFAULTS, settings


class _Default(str):
    """Distinguish an omitted documented default from an explicit setting."""


AUTO = _Default("auto")
STANDARD = _Default("standard")
PREPARATION = {"enrichment", "expand_studies"}


@dataclass(frozen=True)
class PreparationPolicy:
    enrichment: str = "standard"
    expand_studies: bool = True
    linked_metadata: bool | None = None
    peer_provider: bool | None = None

    def allows(self, operation):
        if self.enrichment == "off":
            return False
        if operation == "linked_metadata":
            return self.enrichment == "standard" and self.linked_metadata is not False
        if operation == "peer_provider":
            return self.enrichment == "standard" and self.peer_provider is not False
        return True


def preset(value):
    if not isinstance(value, str) or value not in {"standard", "curators", "off"}:
        raise ValueError("enrichment must be standard, curators or off")
    return str(value)


def _same(left, right):
    if left is right:
        return True
    if isinstance(left, (str, PathLike)) and isinstance(right, (str, PathLike)):
        return str(left) == str(right)
    try:
        return bool(left == right)
    except (TypeError, ValueError):
        return False


def combine(left, right, label):
    result = dict(left)
    for key, value in right.items():
        if key in result and not _same(result[key], value):
            raise ValueError(f"Conflicting {label}: {key}")
        result[key] = value
    return result


def aliases(value, names, label):
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    result = {}
    for key, item in value.items():
        result = combine(result, {names.get(key, key): item}, label)
    return result


def input_settings(value):
    return settings(aliases(value, {"matrix_orientation": "orientation"}, "input options"),
                    set().union(*INPUT.values()) | PREPARATION, "input options")


def output_settings(value, target):
    return settings(aliases(value, {"execution_profile": "profile"}, "output options"),
                    OUTPUT[target], "output options")


def normalize(*, target, in_type, enrichment, options, force_in_type, outfile,
              input_manifest, input_options, output_options, runtime_options):
    simple = aliases(options, {"execution_profile": "profile"}, "options")
    inp = input_settings(input_options)
    out = output_settings(output_options, target)
    runtime = settings(runtime_options, RUNTIME_DEFAULTS, "runtime options")
    paths = {"outfile": outfile, "input_manifest": input_manifest}
    for key, value in simple.items():
        if key in paths:
            if paths[key] is not None and not _same(paths[key], value):
                raise ValueError(f"Conflicting {key}")
            paths[key] = value
        elif key in RUNTIME_DEFAULTS:
            runtime = combine(runtime, {key: value}, "runtime options")
        elif key in PREPARATION or key in set().union(*INPUT.values()):
            # Enrich is historically both an import and an export option.
            if key == "enrich" and key in OUTPUT[target] and key not in inp:
                out = combine(out, {key: value}, "output options")
            else:
                inp = combine(inp, {key: value}, "input options")
        elif key in OUTPUT[target]:
            out = combine(out, {key: value}, "output options")
        else:
            raise ValueError(f"Unsupported options: {key}")
    selected = None if in_type == "auto" else kind_name(in_type)
    old = kind_name(force_in_type)
    if in_type is not AUTO and force_in_type is not None and selected != old:
        raise ValueError("Conflicting in_type and force_in_type")
    if enrichment is not STANDARD:
        inp = combine(inp, {"enrichment": preset(enrichment)}, "preparation options")
    inp.setdefault("enrichment", STANDARD)
    return (selected or old, paths["outfile"], paths["input_manifest"],
            input_settings(inp), output_settings(out, target),
            settings(runtime, RUNTIME_DEFAULTS, "runtime options"))


def preparation_options(inp, out):
    """Extract preparation flags before validating a reader's own settings."""
    inp, out = dict(inp), dict(out)
    configured = inp.pop("enrichment", STANDARD)
    name = preset(configured)
    explicit = configured is not STANDARD
    expand = inp.pop("expand_studies", None)
    old_expand = inp.get("related_series")
    if expand is not None and old_expand is not None and expand != old_expand:
        raise ValueError("Conflicting expand_studies and related_series")
    expand = expand if expand is not None else old_expand if old_expand is not None else True
    if not isinstance(expand, bool):
        raise TypeError("expand_studies must be a boolean")
    legacy = combine({k: v for k, v in inp.items() if k == "enrich"},
                     {k: v for k, v in out.items() if k == "enrich"}, "enrichment options")
    if "enrich" in legacy:
        enabled = legacy["enrich"]
        if explicit and enabled != (name != "off"):
            raise ValueError("Conflicting enrichment and enrich")
        if not explicit:
            name = "standard" if enabled else "off"
    linked, peer = inp.get("enrich_from_geo_ae"), inp.get("include_peer")
    for value in (linked, peer):
        if value is not None and explicit and value != (name == "standard"):
            raise ValueError("Conflicting enrichment preset and repository override")
    return inp, out, PreparationPolicy(name, expand, linked, peer)
