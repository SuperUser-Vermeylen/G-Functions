from g_functions.g_functions import gfunctions


def test_gfunctions_return_shape() -> None:
    result = gfunctions(445)
    assert len(result) == 11
    assert result[0].shape == result[1].shape == result[2].shape == (999,)


def test_phase_diagram_has_required_columns() -> None:
    from g_functions.phase_diagram import build_phase_diagram

    phase_diagram = build_phase_diagram()
    expected = {
        "alphaBeta_L",
        "alphaBeta_R",
        "alphaLiquid_L",
        "alphaLiquid_R",
        "betaLiquid_L",
        "betaLiquid_R",
    }
    assert set(phase_diagram.columns) == expected
    assert len(phase_diagram) > 0


def test_chemical_potential_works_for_valid_temperature() -> None:
    from g_functions.g_functions import chem_potential
    from g_functions.phase_diagram import build_phase_diagram

    phase_diagram = build_phase_diagram()
    g_total_fcc = (1 - (0.1 + 0.2)) * 0
    assert isinstance(chem_potential(400, g_total_fcc, g_total_fcc, g_total_fcc, phase_diagram), list)
