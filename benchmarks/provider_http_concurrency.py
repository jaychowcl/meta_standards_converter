# =============================================================================
# Authors
#
# Created by jaychowcl @ Saez-Rodriguez Group & EMBL-EBI Functional Genomics Team on May 2026
# https://github.com/jaychowcl
# https://saezlab.org
# https://www.ebi.ac.uk/about/teams/functional-genomics/
# =============================================================================
"""Run bounded, model-free concurrency benchmarks against public metadata APIs."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import statistics
import time
from typing import Any, Callable

import requests

from meta_standards_converter.helpers.request_helper import (
    RateLimitedRequester,
    RequestSettings,
)


@dataclass(frozen=True)
class ProviderScenario:
    """One fixed provider endpoint and its conservative client ceiling."""

    name: str
    service: str
    url: str
    items: tuple[dict[str, str], ...]
    request_delay: float
    max_in_flight: int = 2


SCENARIOS = {
    "ncbi_eutils": ProviderScenario(
        name="ncbi_eutils",
        service="ncbi_eutils",
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        items=tuple(
            {"db": "pubmed", "id": pubmed_id, "retmode": "json"}
            for pubmed_id in ("27708338", "31452104", "34216554")
        ),
        request_delay=0.5,
    ),
    "ena_portal": ProviderScenario(
        name="ena_portal",
        service="ena_portal",
        url="https://www.ebi.ac.uk/ena/portal/api/filereport",
        items=tuple(
            {
                "accession": accession,
                "result": "read_run",
                "fields": "run_accession,study_accession",
                "format": "json",
            }
            for accession in ("SRR390728", "SRR341578", "ERR164407")
        ),
        request_delay=1.0,
    ),
    "biostudies": ProviderScenario(
        name="biostudies",
        service="biostudies",
        url="https://www.ebi.ac.uk/biostudies/api/v1/studies/",
        items=tuple(
            {"_path": accession}
            for accession in ("E-MTAB-5061", "E-MTAB-5214", "E-MTAB-6701")
        ),
        request_delay=1.0,
    ),
}


def _percentile_95(values: list[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(0.95 * len(ordered)))]


def _run_mode(
    scenario: ProviderScenario,
    *,
    workers: int,
    requester_factory: Callable[..., Any] = RateLimitedRequester,
) -> dict[str, Any]:
    RateLimitedRequester.reset_service_state()
    settings = RequestSettings(
        request_delay=scenario.request_delay,
        max_in_flight=scenario.max_in_flight,
        max_retries=0,
    )
    starts: list[float] = []

    def fetch(params: dict[str, str]) -> dict[str, Any]:
        request_params = dict(params)
        request_url = scenario.url + request_params.pop("_path", "")

        def measured_get(url: str, **kwargs: Any):
            starts.append(time.monotonic())
            return requests.get(url, **kwargs)

        requester = requester_factory(
            service=scenario.service, settings=settings, get=measured_get
        )
        started = time.monotonic()
        response = requester.get(request_url, params=request_params or None)
        elapsed = time.monotonic() - started
        response.raise_for_status()
        return {
            "status": response.status_code,
            "sha256": hashlib.sha256(response.content).hexdigest(),
            "elapsed_seconds": elapsed,
        }

    mode_started = time.monotonic()
    if workers == 1:
        responses = [fetch(params) for params in scenario.items]
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            responses = list(executor.map(fetch, scenario.items))
    duration = time.monotonic() - mode_started
    ordered_starts = sorted(starts)
    gaps = [right - left for left, right in zip(ordered_starts, ordered_starts[1:])]
    return {
        "duration_seconds": duration,
        "p95_request_seconds": _percentile_95(
            [record["elapsed_seconds"] for record in responses]
        ),
        "minimum_start_gap_seconds": min(gaps) if gaps else None,
        "statuses": [record["status"] for record in responses],
        "digests": [record["sha256"] for record in responses],
    }


def benchmark_scenario(
    scenario: ProviderScenario,
    *,
    runs: int = 3,
    requester_factory: Callable[..., Any] = RateLimitedRequester,
) -> dict[str, Any]:
    """Warm once, compare sequential/two-worker modes, and decide conservatively."""

    _run_mode(scenario, workers=1, requester_factory=requester_factory)
    sequential = [
        _run_mode(scenario, workers=1, requester_factory=requester_factory)
        for _ in range(runs)
    ]
    concurrent = [
        _run_mode(scenario, workers=2, requester_factory=requester_factory)
        for _ in range(runs)
    ]
    sequential_median = statistics.median(
        record["duration_seconds"] for record in sequential
    )
    concurrent_median = statistics.median(
        record["duration_seconds"] for record in concurrent
    )
    speedup = (
        (sequential_median - concurrent_median) / sequential_median
        if sequential_median
        else 0.0
    )
    deterministic = all(
        record["digests"] == sequential[0]["digests"]
        for record in sequential + concurrent
    )
    statuses_ok = all(
        all(status < 400 for status in record["statuses"])
        for record in sequential + concurrent
    )
    sequential_p95 = statistics.median(
        record["p95_request_seconds"] for record in sequential
    )
    concurrent_p95 = statistics.median(
        record["p95_request_seconds"] for record in concurrent
    )
    tail_ok = not sequential_p95 or concurrent_p95 <= sequential_p95 * 1.5
    gap_floor = scenario.request_delay * 0.9
    cap_ok = all(
        record["minimum_start_gap_seconds"] is None
        or record["minimum_start_gap_seconds"] >= gap_floor
        for record in concurrent
    )
    enabled = speedup >= 0.20 and deterministic and statuses_ok and tail_ok and cap_ok
    return {
        "scenario": asdict(scenario),
        "method": {"warmups": 1, "runs": runs, "workers": [1, 2]},
        "sequential": sequential,
        "concurrent": concurrent,
        "decision": {
            "enable_two_workers": enabled,
            "speedup_fraction": speedup,
            "deterministic": deterministic,
            "statuses_ok": statuses_ok,
            "tail_latency_ok": tail_ok,
            "start_rate_cap_ok": cap_ok,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--provider", choices=(*SCENARIOS, "all"), default="all"
    )
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.runs < 1:
        parser.error("--runs must be positive")
    selected = SCENARIOS.values() if args.provider == "all" else (SCENARIOS[args.provider],)
    report = {
        "schema_version": "1.0",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "model_calls": 0,
        "embedding_calls": 0,
        "results": [benchmark_scenario(item, runs=args.runs) for item in selected],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
