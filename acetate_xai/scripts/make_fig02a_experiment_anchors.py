from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Re-generate Fig02A experimental anchors (4-panel) from conditions CSV.")
    p.add_argument(
        "--conditions-csv",
        default=str(Path("acetate_xai") / "data" / "conditions_experiment.csv"),
        help="Path to conditions_experiment.csv (must include measured_OD).",
    )
    p.add_argument(
        "--outdir",
        default=str(Path("results") / "figures_draft"),
        help="Output directory (default: results/figures_draft).",
    )
    p.add_argument("--figsize-x", type=int, default=10, help="Figure width in inches (default: 10).")
    p.add_argument("--figsize-y", type=int, default=7, help="Figure height in inches (default: 7).")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Reuse the exact plotting implementation used in the paper figure pipeline.
    # NOTE: This import is intentionally local to avoid importing heavy modules when just showing --help.
    from make_paper_figures import _plot_experiment_anchors_4panel  # noqa: PLC0415

    _plot_experiment_anchors_4panel(
        conditions_csv=Path(args.conditions_csv),
        out_png=outdir / "Fig02A_experiment_anchors_4panel.png",
        out_csv=outdir / "Fig02A_experiment_anchors_4panel.csv",
        figsize=(int(args.figsize_x), int(args.figsize_y)),
    )
    print(f"[OK] Wrote: {outdir / 'Fig02A_experiment_anchors_4panel.png'}")
    print(f"[OK] Wrote: {outdir / 'Fig02A_experiment_anchors_4panel.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

