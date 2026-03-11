# Testing calculation scheme generation through case_generator
# R.Guasch
# Date : 10/03/2026

import pytest
from starterDD.InterfaceToDD.case_generator import DragonCase
from starterDD.InterfaceToDD.CLE2000 import (
    main_procedure, sub_procedure,
    CLE2000_MAX_LINE, CLE2000_MAX_VARNAME,
    validate_varname, wrap_cle2000_line,
)
from starterDD.DDModel.DragonCalculationScheme import (
    CalculationStep, DragonCalculationScheme,
)


# ---------------------------------------------------------
# Test CLE2000 constraints
# ---------------------------------------------------------

def test_varname_length_limit():
    """Variable names > 12 chars must raise."""
    with pytest.raises(ValueError, match="12"):
        validate_varname("this_is_too_long")
    # Exactly 12 is fine
    validate_varname("twelve_chars")


def test_line_wrapping():
    """Lines > 70 chars get wrapped."""
    long = "A " * 40  # 80 chars
    wrapped = wrap_cle2000_line(long.strip())
    for line in wrapped.splitlines():
        assert len(line) <= CLE2000_MAX_LINE


# ---------------------------------------------------------
# Integration: full case generation
# ---------------------------------------------------------

def test_generate_cle2000_procedures(tmp_path):
    """
    Generate full CLE2000 case from GE14 DOM config yamls
    and verify all four files are produced.
    """
    GE14_DOM_test_case = DragonCase(
        case_name="GE14_DOM",
        call_glow=False,
        draglibs_names_to_alias={
            "draglibendfb8r1SHEM295": "endfb8r1_295",
        },
        config_yamls={
            "MATS": "data/BWRProgressionProblems"
                    "/GE14/inputs"
                    "/material_compositions.yaml",
            "GEOM": "data/BWRProgressionProblems"
                    "/GE14/inputs"
                    "/GEOM_GE14_DOM.yaml",
            "CALC_SCHEME": "data/BWRProgressionProblems"
                           "/GE14/inputs"
                           "/CALC_SCHEME_1L.yaml",
        },
        output_path=str(tmp_path),
    )

    result = GE14_DOM_test_case.generate_cle2000_procedures()

    # All four files should exist
    assert os.path.isfile(result["x2m"])
    assert os.path.isfile(result["mix"])
    assert os.path.isfile(result["trk"])
    assert os.path.isfile(result["edir"])

    # Verify CLE2000 line-length constraint on x2m
    with open(result["x2m"]) as f:
        for i, line in enumerate(f, 1):
            stripped = line.rstrip("\n")
            assert len(stripped) <= CLE2000_MAX_LINE, (
                f"x2m line {i} exceeds {CLE2000_MAX_LINE} "
                f"chars: {stripped!r}"
            )


import os