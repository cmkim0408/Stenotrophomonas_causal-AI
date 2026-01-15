from __future__ import annotations

from typer.testing import CliRunner

from acetate_xai.cli import app
from acetate_xai.io import REQUIRED_CONDITIONS_COLUMNS, load_conditions_csv
from acetate_xai.medium import apply_condition_to_model
from acetate_xai.config import load_config


def test_cli_help_runs() -> None:
    runner = CliRunner()
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "acetate_xai" in result.stdout


def test_load_conditions_example_csv() -> None:
    df = load_conditions_csv("data/conditions_experiment_example.csv")
    for c in REQUIRED_CONDITIONS_COLUMNS:
        assert c in df.columns
    assert len(df) >= 5


def test_apply_condition_to_model_changes_bounds() -> None:
    """
    Minimal unit test without a full GSM:
    create a toy model with 2 exchange reactions and verify bounds change.
    """
    from cobra import Model, Reaction, Metabolite

    m = Model("toy")
    met = Metabolite("a_e", compartment="e")

    ex_ac = Reaction("EX_ac_e")
    ex_ac.name = "acetate exchange"
    ex_ac.add_metabolites({met: -1})
    ex_ac.lower_bound = 0.0
    ex_ac.upper_bound = 1000.0

    ex_nh4 = Reaction("EX_nh4_e")
    ex_nh4.name = "ammonium exchange"
    ex_nh4.add_metabolites({met: -1})
    ex_nh4.lower_bound = 0.0
    ex_nh4.upper_bound = 1000.0

    m.add_reactions([ex_ac, ex_nh4])

    cfg = load_config("configs/medium_base.yaml")
    row = {
        "condition_id": "toy",
        "pH0": 7.0,
        "acetate_mM": 100.0,
        "nh4cl_gL": 1.0,
        "yeast_extract_gL": 0.0,
    }
    res = apply_condition_to_model(m, row, cfg)
    assert res.condition_id == "toy"
    assert m.reactions.get_by_id("EX_ac_e").lower_bound < 0
    assert m.reactions.get_by_id("EX_nh4_e").lower_bound < 0

