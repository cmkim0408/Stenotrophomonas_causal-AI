from __future__ import annotations

import logging

import typer
from rich.logging import RichHandler

from acetate_xai import __version__

app = typer.Typer(add_completion=False, help="acetate_xai: GSM simulation + XAI pipeline")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, show_time=False, show_path=False)],
    )


@app.callback()
def main(
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable verbose logging."),
) -> None:
    """Entry point."""
    _setup_logging(verbose=verbose)


@app.command()
def version() -> None:
    """Print package version."""
    typer.echo(__version__)


if __name__ == "__main__":
    app()

