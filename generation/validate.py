"""Candidate fixture validation and pipeline consumer inventory."""

from __future__ import annotations

import csv
import json
import re
import subprocess
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .core import STAGES, h5ad_compression_violations, sha256_file

TEST_PATTERN = re.compile(r"\btest\(\s*(['\"])(.*?)\1", re.DOTALL)
PIPELINE_ROUTE = re.compile(
    r"(?:(?:params\.)?pipelines_testdata_base_path|fixtureBase)\s*\+\s*(['\"])(.*?)\1",
    re.DOTALL,
)
ROOTED_TEXT_ROUTE = re.compile(r"\$\{fixtureBase\}([^,\s'\"\\]+)")
MODULE_ROUTE = re.compile(r"params\.modules_testdata_base_path\s*\+\s*(['\"])(.*?)\1", re.DOTALL)
URL_ROUTE = re.compile(r"https?://[^\s'\"`),]+")


@dataclass(frozen=True)
class ConsumerCase:
    """One CI-enabled nf-test case and the fixture routes visible in its source."""

    path: str
    name: str
    category: str | None
    behaviour: str
    expected_failure: bool
    route_classifications: list[str]
    routes: list[str]
    unresolved: list[str]


def audit_consumers(
    pipeline_dir: Path, contracts: dict[str, Any], dataset: Path | None = None
) -> dict[str, Any]:
    """Audit every local CI case and the CSV routes that can bypass fixture routing."""
    pipeline_dir = pipeline_dir.resolve()
    audit_contract = contracts["consumer_audit"]
    required_categories = set(audit_contract["required_categories"])
    dynamic_patterns = set(audit_contract["dynamic_dataset_patterns"])
    cases: list[ConsumerCase] = []
    external_fixture_routes: set[str] = set()

    for root in audit_contract["included_roots"]:
        directory = pipeline_dir / root
        if not directory.exists():
            raise ValueError(f"Audit root does not exist: {directory}")
        for test_path in sorted(directory.rglob("*.nf.test")):
            text = test_path.read_text(encoding="utf-8")
            external_fixture_routes.update(
                match.group(2).lstrip("/") for match in MODULE_ROUTE.finditer(text)
            )
            category = _suite_category(text, required_categories)
            first_case = TEST_PATTERN.search(text)
            shared_setup = text[: first_case.start()] if first_case else ""
            for name, body in _test_blocks(text):
                if category is None:
                    continue
                cases.append(
                    _classify_case(
                        test_path.relative_to(pipeline_dir).as_posix(),
                        name,
                        category,
                        shared_setup + body,
                        dynamic_patterns,
                    )
                )

    csv_routes = _audit_csv_routes(pipeline_dir)
    unresolved = [
        {"path": case.path, "name": case.name, "expressions": case.unresolved}
        for case in cases
        if case.unresolved
    ]
    expected_cases = audit_contract["expected_cases"]
    category_counts = Counter(case.category for case in cases)
    behaviour_counts = Counter(case.behaviour for case in cases)
    route_counts = Counter(
        route_classification
        for case in cases
        for route_classification in case.route_classifications
    )
    missing_candidate_routes: list[dict[str, str]] = []
    if dataset is not None:
        dataset = dataset.resolve()
        for case in cases:
            for route in case.routes:
                if route.startswith(("http://", "https://")) or route in external_fixture_routes:
                    continue
                pattern = re.sub(r"\$\{[^}]+\}", "*", route)
                matches = [path for path in dataset.glob(pattern) if path.is_file()]
                if not matches:
                    missing_candidate_routes.append(
                        {
                            "suite": case.path,
                            "test": case.name,
                            "route": route,
                        }
                    )
    report = {
        "schema_version": 1,
        "pipeline": str(pipeline_dir),
        "pipeline_git_state": _git_state(pipeline_dir),
        "expected_case_count": expected_cases,
        "historical_baseline_case_count": audit_contract.get("baseline_cases"),
        "baseline_reconciliation": audit_contract.get("baseline_reconciliation"),
        "actual_case_count": len(cases),
        "case_count_matches_contract": len(cases) == expected_cases,
        "category_counts": dict(sorted(category_counts.items())),
        "behaviour_counts": dict(sorted(behaviour_counts.items())),
        "route_classification_counts": dict(sorted(route_counts.items())),
        "cases": [asdict(case) for case in cases],
        "csv_routes": csv_routes,
        "unresolved_classifications": unresolved,
        "candidate": str(dataset) if dataset is not None else None,
        "missing_candidate_routes": missing_candidate_routes,
        "coverage_status": "complete"
        if len(cases) == expected_cases and not unresolved and not missing_candidate_routes
        else "incomplete",
        "execution_evidence": (
            "This is a static consumer inventory. Full nf-test category runs remain required "
            "to establish runtime coverage."
        ),
    }
    return report


def validate_dataset(
    dataset: Path,
    contracts: dict[str, Any],
    family: str | None = None,
    pipeline_revision: str | None = None,
) -> dict[str, Any]:
    """Validate size and H5AD compression gates for a candidate fixture tree."""
    dataset = dataset.resolve()
    selected_families = _selected_families(contracts, family)
    files: list[Path] = []
    missing_patterns: list[str] = []
    declared_patterns: list[str] = []
    for fixture_family in selected_families:
        for pattern in contracts["fixture_families"][fixture_family]["paths"]:
            declared_patterns.append(pattern)
            matches = sorted(dataset.glob(pattern))
            if not matches:
                missing_patterns.append(pattern)
            files.extend(match for match in matches if match.is_file())

    size_policy = contracts["size_policy"]
    hard_limit = size_policy["hard_limit_bytes"]
    target = size_policy["target_bytes"]
    hard_limit_violations: list[dict[str, Any]] = []
    target_reviews: list[dict[str, Any]] = []
    compression_violations: list[dict[str, Any]] = []
    if family is None:
        files.extend(path for path in dataset.rglob("*") if path.is_file())
    seen: set[Path] = set()
    for path in files:
        if path in seen:
            continue
        seen.add(path)
        size = path.stat().st_size
        relative = path.relative_to(dataset).as_posix()
        if size >= hard_limit:
            hard_limit_violations.append({"path": relative, "bytes": size})
        elif size >= target:
            target_reviews.append({"path": relative, "bytes": size})
        if path.suffix == ".h5ad":
            invalid_datasets = h5ad_compression_violations(path)
            if invalid_datasets:
                compression_violations.append({"path": relative, "datasets": invalid_datasets})

    expected_revision = contracts["pipeline"]["revision"]
    revision_matches_contract = pipeline_revision in {None, expected_revision}
    undeclared_files = []
    if family is None:
        undeclared_files = sorted(
            path.relative_to(dataset).as_posix()
            for path in seen
            if not _matches_any(path.relative_to(dataset), declared_patterns)
            and not path.relative_to(dataset).as_posix().startswith("reports/")
        )
    report_errors: list[str] = []
    if family is None:
        for stage in STAGES:
            report_path = dataset / "reports" / f"{stage}.json"
            if not report_path.is_file():
                report_errors.append(f"Missing stage report: {stage}")
                continue
            try:
                stage_report = json.loads(report_path.read_text(encoding="utf-8"))
            except (ValueError, OSError) as error:
                report_errors.append(f"Invalid stage report {stage}: {error}")
                continue
            produced = stage_report.get("produced", [])
            if not produced or len(set(produced)) != len(produced):
                report_errors.append(f"Stage report {stage} has no unique output inventory")
            for relative in produced:
                path = (dataset / relative).resolve()
                if not path.is_relative_to(dataset) or not path.is_file():
                    report_errors.append(f"Invalid reported output: {relative}")
                    continue
                record = stage_report.get("artefacts", {}).get(relative, {})
                if record.get("sha256") != sha256_file(path):
                    report_errors.append(f"Reported output checksum mismatch: {relative}")
                if record.get("bytes") != path.stat().st_size:
                    report_errors.append(f"Reported output size mismatch: {relative}")
    passed = (
        not report_errors
        and not missing_patterns
        and not hard_limit_violations
        and not compression_violations
        and revision_matches_contract
        and not undeclared_files
    )
    return {
        "schema_version": 1,
        "dataset": str(dataset),
        "families": selected_families,
        "files_checked": len(seen),
        "missing_patterns": missing_patterns,
        "undeclared_files": undeclared_files,
        "hard_limit_violations": hard_limit_violations,
        "target_size_reviews": target_reviews,
        "h5ad_compression_violations": compression_violations,
        "provenance": {
            "expected_pipeline_revision": expected_revision,
            "actual_pipeline_revision": pipeline_revision,
            "pipeline_revision_matches_contract": revision_matches_contract,
        },
        "stage_report_errors": report_errors,
        "status": "passed" if passed else "failed",
    }


def write_report(report: dict[str, Any], path: Path) -> None:
    """Write a deterministic JSON report for review and later publication."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _test_blocks(text: str) -> list[tuple[str, str]]:
    blocks: list[tuple[str, str]] = []
    for match in TEST_PATTERN.finditer(text):
        opening = text.find("{", match.end())
        if opening == -1:
            continue
        closing = _matching_brace(text, opening)
        if closing is None:
            continue
        blocks.append((match.group(2), text[opening + 1 : closing]))
    return blocks


def _matching_brace(text: str, opening: int) -> int | None:
    depth = 0
    quote: str | None = None
    escaped = False
    index = opening
    while index < len(text):
        character = text[index]
        if quote is not None:
            if len(quote) == 3 and text.startswith(quote, index):
                index += 3
                quote = None
                continue
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == quote:
                quote = None
            index += 1
            continue
        if text.startswith("//", index):
            newline = text.find("\n", index + 2)
            index = len(text) if newline < 0 else newline + 1
            continue
        if text.startswith("/*", index):
            closing = text.find("*/", index + 2)
            if closing < 0:
                return None
            index = closing + 2
            continue
        if text.startswith("'''", index) or text.startswith('"""', index):
            quote = text[index : index + 3]
            index += 3
            continue
        if character in {"'", '"'}:
            quote = character
        elif character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth == 0:
                return index
        index += 1
    return None


def _suite_category(text: str, allowed: set[str]) -> str | None:
    for category in sorted(allowed):
        if re.search(rf"\btag\s+['\"]{re.escape(category)}['\"]", text):
            return category
    return None


def _classify_case(
    path: str, name: str, category: str, body: str, dynamic_patterns: set[str]
) -> ConsumerCase:
    routes: list[str] = []
    route_classifications: list[str] = []
    unresolved: list[str] = []
    for match in PIPELINE_ROUTE.finditer(body):
        route = match.group(2).replace("\\$", "$").lstrip("/")
        routes.append(route)
        route_classifications.append("dataset")
        if "${" in route and route not in dynamic_patterns:
            unresolved.append(route)
    for match in ROOTED_TEXT_ROUTE.finditer(body):
        route = match.group(1).lstrip("/")
        routes.append(route)
        route_classifications.append("dataset")
    for match in MODULE_ROUTE.finditer(body):
        routes.append(match.group(2).lstrip("/"))
        route_classifications.append("modules_external")
    for route in URL_ROUTE.findall(body):
        routes.append(route)
        route_classifications.append("hardcoded_external")

    if not routes:
        route_classifications.append("inline_or_derived")
    lowered = name.lower()
    expected_failure = bool(re.search(r"\b(?:process|workflow)\.failed\b", body))
    behaviour = _behaviour_class(lowered, expected_failure)
    return ConsumerCase(
        path=path,
        name=name,
        category=category,
        behaviour=behaviour,
        expected_failure=expected_failure,
        route_classifications=sorted(set(route_classifications)),
        routes=sorted(set(routes)),
        unresolved=sorted(set(unresolved)),
    )


def _behaviour_class(name: str, expected_failure: bool) -> str:
    if expected_failure or "fail" in name or "reject" in name:
        return "expected_failure"
    structural_terms = ("duplicate", "isoform", "column", "layer", "matrix", "missing")
    if any(term in name for term in structural_terms):
        return "structural"
    if any(term in name for term in ("plot", "snapshot", "report")):
        return "reporting_or_plotting"
    numerical_terms = ("mapping", "extension", "integrat", "differential", "normalis")
    if any(term in name for term in numerical_terms):
        return "biological_or_numerical"
    return "behavioural"


def _audit_csv_routes(pipeline_dir: Path) -> list[dict[str, Any]]:
    routes: list[dict[str, Any]] = []
    for csv_path in sorted((pipeline_dir / "tests").glob("*.csv")):
        with csv_path.open(encoding="utf-8", newline="") as handle:
            for row_number, row in enumerate(csv.DictReader(handle), start=2):
                for column, value in row.items():
                    if value and URL_ROUTE.fullmatch(value):
                        routes.append(
                            {
                                "path": csv_path.relative_to(pipeline_dir).as_posix(),
                                "row": row_number,
                                "column": column,
                                "value": value,
                                "classification": "hardcoded_external",
                            }
                        )
    return routes


def _selected_families(contracts: dict[str, Any], family: str | None) -> list[str]:
    families = contracts["fixture_families"]
    if family is None:
        return sorted(families)
    if family in families:
        return [family]
    if family in STAGES:
        return sorted(name for name, spec in families.items() if spec["stage"] == family)
    raise ValueError(f"Unknown family or stage: {family}")


def _matches_any(path: Path, patterns: list[str]) -> bool:
    return any(path.match(pattern) for pattern in patterns)


def _git_state(directory: Path) -> dict[str, Any]:
    """Record whether the inventory came from an immutable pipeline checkout."""
    revision = subprocess.check_output(
        ["git", "-C", str(directory), "rev-parse", "HEAD"], text=True
    ).strip()
    changed_paths = subprocess.check_output(
        ["git", "-C", str(directory), "status", "--short"], text=True
    ).splitlines()
    return {
        "revision": revision,
        "working_tree_dirty": bool(changed_paths),
        "changed_paths": changed_paths,
    }
