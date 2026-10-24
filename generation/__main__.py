"""Command-line interface for building, validating and auditing fixtures."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable
from pathlib import Path
from time import perf_counter

from . import analysis, reference, singler, source, structural
from .core import (
    STAGES,
    ContractError,
    StageNotImplementedError,
    assert_stage_output_is_available,
    create_context,
    load_contracts,
    pipeline_revision,
    write_stage_report,
)
from .validate import audit_consumers, validate_dataset, write_report

BUILDERS: dict[str, Callable] = {
    "source": source.build,
    "structural": structural.build,
    "analysis": analysis.build,
    "reference": reference.build,
    "singler": singler.build,
}


def main(argv: list[str] | None = None) -> int:
    """Run the fixture-builder command line."""
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            return _build(args)
        if args.command == "validate":
            return _validate(args)
        if args.command == "audit":
            return _audit(args)
    except (ContractError, StageNotImplementedError, OSError, ValueError) as error:
        parser.error(str(error))
    return 2


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)

    build = subcommands.add_parser("build", help="build one declared stage")
    build.add_argument("--stage", choices=[*STAGES, "all"], required=True)
    build.add_argument("--out", type=Path, required=True, help="external candidate directory")
    build.add_argument("--cache", type=Path, required=True, help="source download cache")
    build.add_argument("--pipeline", type=Path, required=True)
    build.add_argument("--contracts", type=Path, default=Path("fixtures.json"))

    validate = subcommands.add_parser("validate", help="validate a candidate fixture tree")
    validate.add_argument("--dataset", type=Path, required=True)
    validate.add_argument("--pipeline", type=Path, required=True)
    validate.add_argument("--report", type=Path, required=True)
    validate.add_argument("--contracts", type=Path, default=Path("fixtures.json"))
    validate.add_argument("--family")

    audit = subcommands.add_parser("audit", help="inventory all local CI fixture consumers")
    audit.add_argument("--pipeline", type=Path, required=True)
    audit.add_argument("--contracts", type=Path, default=Path("fixtures.json"))
    audit.add_argument("--report", type=Path)
    audit.add_argument("--dataset", type=Path, help="also verify candidate fixture routes")
    return parser


def _build(args: argparse.Namespace) -> int:
    stages = STAGES if args.stage == "all" else (args.stage,)
    reports: list[str] = []
    for stage in stages:
        context = create_context(
            stage=stage,
            contracts_path=args.contracts,
            out_dir=args.out,
            cache_dir=args.cache,
            pipeline_dir=args.pipeline,
        )
        assert_stage_output_is_available(context)
        started = perf_counter()
        result = BUILDERS[stage](context)
        reports.append(str(write_stage_report(context, result, started)))
    print(json.dumps({"status": "passed", "reports": reports}, sort_keys=True))
    return 0


def _validate(args: argparse.Namespace) -> int:
    contracts = load_contracts(args.contracts.resolve())
    resolved_pipeline = args.pipeline.resolve()
    report = validate_dataset(
        args.dataset,
        contracts,
        args.family,
        pipeline_revision=pipeline_revision(resolved_pipeline),
    )
    report["pipeline"] = str(resolved_pipeline)
    write_report(report, args.report)
    print(json.dumps({"report": str(args.report), "status": report["status"]}, sort_keys=True))
    return 0 if report["status"] == "passed" else 1


def _audit(args: argparse.Namespace) -> int:
    contracts = load_contracts(args.contracts.resolve())
    report = audit_consumers(args.pipeline, contracts, dataset=args.dataset)
    if args.report is not None:
        write_report(report, args.report)
    print(
        json.dumps(
            {
                "cases": report["actual_case_count"],
                "coverage_status": report["coverage_status"],
                "report": str(args.report) if args.report else None,
                "unresolved": len(report["unresolved_classifications"]),
            },
            sort_keys=True,
        )
    )
    return 0 if report["coverage_status"] == "complete" else 1


if __name__ == "__main__":
    sys.exit(main())
