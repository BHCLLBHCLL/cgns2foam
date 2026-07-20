"""Templates for chtMultiRegionSimpleFoam case files.

Holds the per-region thermophysicalProperties / turbulenceProperties /
radiationProperties, fvSchemes / fvSolution / fvOptions, 0/ field stubs,
controlDict, MRFProperties and decomposeParDict builders used by
:mod:`src.cht_direct` (the one-step multi-region writer).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .regions_config import MaterialSpec


# ---------------------------------------------------------------------------
# Minimal OpenFOAM dictionary helpers
# ---------------------------------------------------------------------------


def _foam_header(obj_class: str, obj_name: str, location: str = "") -> str:
    loc = f'\n    location    "{location}";' if location else ""
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2412                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {obj_class};{loc}
    object      {obj_name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8", newline="\n")


# ---------------------------------------------------------------------------
# Default material / numerics
# ---------------------------------------------------------------------------

_DEFAULT_FLUID_MAT: dict[str, Any] = {
    "thermoType": {
        "type": "heRhoThermo",
        "mixture": "pureMixture",
        "transport": "const",
        "thermo": "hConst",
        "equationOfState": "perfectGas",
        "specie": "specie",
        "energy": "sensibleEnthalpy",
    },
    "mixture": {
        "specie": {"nMoles": 1, "molWeight": 28.966},
        "thermodynamics": {"Cp": 1006.43, "Hf": 0},
        "transport": {"mu": 1.846e-5, "Pr": 0.706},
    },
}

_DEFAULT_SOLID_MAT: dict[str, Any] = {
    "thermoType": {
        "type": "heSolidThermo",
        "mixture": "pureMixture",
        "transport": "constIso",
        "thermo": "hConst",
        "equationOfState": "rhoConst",
        "specie": "specie",
        "energy": "sensibleEnthalpy",
    },
    "mixture": {
        "specie": {"nMoles": 1, "molWeight": 26.98},
        "thermodynamics": {"Hf": 0, "Sf": 0, "Cp": 871},
        "transport": {"kappa": 202.4},
        "equationOfState": {"rho": 2719},
    },
}


def _thermophysical_fluid(mat: MaterialSpec | None = None) -> str:
    t = _DEFAULT_FLUID_MAT["thermoType"]
    m = {
        "specie": dict(_DEFAULT_FLUID_MAT["mixture"]["specie"]),
        "thermodynamics": dict(_DEFAULT_FLUID_MAT["mixture"]["thermodynamics"]),
        "transport": dict(_DEFAULT_FLUID_MAT["mixture"]["transport"]),
    }
    if mat is not None:
        if mat.mol_weight is not None:
            m["specie"]["molWeight"] = mat.mol_weight
        if mat.cp is not None:
            m["thermodynamics"]["Cp"] = mat.cp
        if mat.mu is not None:
            m["transport"]["mu"] = mat.mu
        if mat.pr is not None:
            m["transport"]["Pr"] = mat.pr
    return (
        _foam_header("dictionary", "thermophysicalProperties", "constant")
        + f"""
thermoType
{{
    type            {t['type']};
    mixture         {t['mixture']};
    transport       {t['transport']};
    thermo          {t['thermo']};
    equationOfState {t['equationOfState']};
    specie          {t['specie']};
    energy          {t['energy']};
}}

mixture
{{
    specie
    {{
        nMoles          {m['specie']['nMoles']};
        molWeight       {m['specie']['molWeight']};
    }}
    thermodynamics
    {{
        Cp              {m['thermodynamics']['Cp']};
        Hf              {m['thermodynamics']['Hf']};
    }}
    transport
    {{
        mu              {m['transport']['mu']};
        Pr              {m['transport']['Pr']};
    }}
}}

// ************************************************************************* //
"""
    )


def _thermophysical_solid(mat: MaterialSpec | None = None) -> str:
    t = _DEFAULT_SOLID_MAT["thermoType"]
    m = {
        "specie": dict(_DEFAULT_SOLID_MAT["mixture"]["specie"]),
        "thermodynamics": dict(_DEFAULT_SOLID_MAT["mixture"]["thermodynamics"]),
        "transport": dict(_DEFAULT_SOLID_MAT["mixture"]["transport"]),
        "equationOfState": dict(_DEFAULT_SOLID_MAT["mixture"]["equationOfState"]),
    }
    if mat is not None:
        if mat.mol_weight is not None:
            m["specie"]["molWeight"] = mat.mol_weight
        if mat.cp is not None:
            m["thermodynamics"]["Cp"] = mat.cp
        if mat.kappa is not None:
            m["transport"]["kappa"] = mat.kappa
        if mat.rho is not None:
            m["equationOfState"]["rho"] = mat.rho
    return (
        _foam_header("dictionary", "thermophysicalProperties", "constant")
        + f"""
thermoType
{{
    type            {t['type']};
    mixture         {t['mixture']};
    transport       {t['transport']};
    thermo          {t['thermo']};
    equationOfState {t['equationOfState']};
    specie          {t['specie']};
    energy          {t['energy']};
}}

mixture
{{
    specie
    {{
        nMoles          {m['specie']['nMoles']};
        molWeight       {m['specie']['molWeight']};
    }}
    thermodynamics
    {{
        Hf              {m['thermodynamics']['Hf']};
        Sf              {m['thermodynamics']['Sf']};
        Cp              {m['thermodynamics']['Cp']};
    }}
    transport
    {{
        kappa           {m['transport']['kappa']};
    }}
    equationOfState
    {{
        rho             {m['equationOfState']['rho']};
    }}
}}

// ************************************************************************* //
"""
    )


def _turbulence_laminar() -> str:
    return (
        _foam_header("dictionary", "turbulenceProperties", "constant")
        + """
simulationType laminar;

// ************************************************************************* //
"""
    )


def _control_dict_cht(
    *,
    end_time: int = 500,
    write_interval: int = 50,
    purge_write: int = 0,
) -> str:
    return (
        _foam_header("dictionary", "controlDict")
        + f"""
application     chtMultiRegionSimpleFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end_time};

deltaT          1;

writeControl    timeStep;
writeInterval   {write_interval};
purgeWrite      {purge_write};

writeFormat     binary;
writePrecision  8;
writeCompression off;

timeFormat      general;
timePrecision   8;
runTimeModifiable true;

// Convergence is driven by residualControl in system/<region>/fvSolution;
// per-field residuals are also printed to the solver log.  To monitor
// fields, uncomment e.g.:
//
// functions
// {{
//     // min/max per region (region name as in constant/regionProperties)
//     // #includeFunc fieldMinMax(region=air,fields=(T U p_rgh))
// }}

// ************************************************************************* //
"""
    )


def _fv_schemes_fluid() -> str:
    return (
        _foam_header("dictionary", "fvSchemes", "system")
        + """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }

// First-order upwind for robust start-up.  Once the solution has settled
// (~100-200 iterations), switch div(phi,...) to second order for accuracy:
//
//     div(phi,U)      bounded Gauss linearUpwind grad(U);
//     div(phi,h)      bounded Gauss linearUpwind grad(h);
//     div(phi,K)      bounded Gauss linearUpwind grad(K);
//
divSchemes
{
    div(phi,U)      bounded Gauss upwind;
    div(phi,h)      bounded Gauss upwind;
    div(phi,K)      bounded Gauss upwind;
    div((muEff*dev2(T(grad(U))))) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited 0.333; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.333; }
wallDist { method meshWave; }

// ************************************************************************* //
"""
    )


def _fv_schemes_solid() -> str:
    return (
        _foam_header("dictionary", "fvSchemes", "system")
        + """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default Gauss linear; }
laplacianSchemes { default Gauss linear limited 0.333; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.333; }

// ************************************************************************* //
"""
    )


def _fv_solution_fluid() -> str:
    return (
        _foam_header("dictionary", "fvSolution", "system")
        + """
solvers
{
    p_rgh
    {
        solver           GAMG;
        smoother         GaussSeidel;
        tolerance        1e-7;
        relTol           0.01;
    }
    "(U|h)"
    {
        solver           PBiCGStab;
        preconditioner   DILU;
        tolerance        1e-6;
        relTol           0.05;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 1;
    momentumPredictor true;
    residualControl
    {
        p_rgh           1e-5;
        U               1e-5;
        h               1e-6;
    }
}
relaxationFactors
{
    fields { p_rgh 0.3; rho 0.05; }
    equations { U 0.3; h 0.3; }
}

// ************************************************************************* //
"""
    )


def _fv_solution_solid() -> str:
    return (
        _foam_header("dictionary", "fvSolution", "system")
        + """
solvers
{
    h
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-8;
        relTol          0.05;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 1;
    residualControl
    {
        h               1e-6;
    }
}
relaxationFactors
{
    equations { h 1; }
}

// ************************************************************************* //
"""
    )


def _region_properties(fluid: list[str], solid: list[str]) -> str:
    return (
        _foam_header("dictionary", "regionProperties")
        + f"""
regions
(
    fluid ( {' '.join(fluid)} )
    solid ( {' '.join(solid)} )
);

// ************************************************************************* //
"""
    )


def _fv_options_fluid() -> str:
    return (
        _foam_header("dictionary", "fvOptions", "system")
        + """
limitT
{
    type            limitTemperature;
    active          yes;
    selectionMode   all;
    min             200;
    max             500;
}

limitU
{
    type            limitVelocity;
    active          yes;
    selectionMode   all;
    max             100;
}

// ************************************************************************* //
"""
    )


def _fv_options_solid_heat(power_watts: float) -> str:
    """Volumetric heat source for a solid region (total power in watts).

    ``volumeMode absolute`` makes OpenFOAM interpret ``explicit`` as the
    *total* source integrated over the selected cells (watts for the
    enthalpy equation); the per-volume density is derived from the cell
    volumes internally.
    """
    return (
        _foam_header("dictionary", "fvOptions", "system")
        + f"""
heatSource
{{
    type            scalarSemiImplicitSource;
    active          yes;
    selectionMode   all;
    volumeMode      absolute;
    sources
    {{
        h
        {{
            explicit       {power_watts};
            implicit       0;
        }}
    }}
}}

// ************************************************************************* //
"""
    )


def _radiation_none() -> str:
    return (
        _foam_header("dictionary", "radiationProperties", "constant")
        + """
radiation       off;
radiationModel  none;

// ************************************************************************* //
"""
    )


def _gravity(g: list[float] | None = None) -> str:
    gx, gy, gz = g or [0.0, 0.0, -9.81]
    return (
        _foam_header("uniformDimensionedVectorField", "g", "constant")
        + f"""
dimensions      [0 1 -2 0 0 0 0];
value           ({gx} {gy} {gz});
// ************************************************************************* //
"""
    )


def mrf_non_rotating_patches(patch_names: list[str]) -> list[str]:
    """Patches that stay in the inertial frame under MRF."""
    out: list[str] = []
    for p in patch_names:
        low = p.lower()
        if low.startswith("ami_") or ("ami" in low and "_to_" not in low):
            out.append(p)
        elif p == "open" or p.startswith("open"):
            out.append(p)
        elif "_to_" in p:
            out.append(p)
    return sorted(set(out))


def mrf_properties(
    entries: list[dict[str, Any]],
    *,
    non_rotating: list[str] | None = None,
    location: str = "constant",
) -> str:
    """Build ``MRFProperties`` from resolved MRF entries.

    Each *entries* item needs keys: ``name``, ``cellZone``, ``origin``,
    ``axis``, ``omega``.  Optional per-entry ``nonRotatingPatches`` overrides
    the shared *non_rotating* list.
    """
    shared_nr = non_rotating or []
    blocks: list[str] = []
    for i, e in enumerate(entries):
        name = str(e.get("name") or (f"MRF{i + 1}" if len(entries) > 1 else "MRF"))
        cz = e["cellZone"]
        ox, oy, oz = e["origin"]
        ax, ay, az = e["axis"]
        omega = float(e["omega"])
        nr = e.get("nonRotatingPatches")
        if nr is None:
            nr = shared_nr
        if nr:
            nr_block = f"nonRotatingPatches ( {' '.join(nr)} );"
        else:
            nr_block = "nonRotatingPatches ();"
        blocks.append(
            f"""{name}
{{
    cellZone            {cz};
    active              yes;
    {nr_block}
    origin              ({ox} {oy} {oz});
    axis                ({ax} {ay} {az});
    omega               {omega};
}}"""
        )
    body = "\n\n".join(blocks) if blocks else "// (no MRF zones)"
    return (
        _foam_header("dictionary", "MRFProperties", location)
        + f"""

{body}

// ************************************************************************* //
"""
    )


def _is_cyclic_ami_patch(name: str, patch_types: dict[str, str] | None) -> bool:
    if patch_types and patch_types.get(name) == "cyclicAMI":
        return True
    low = name.lower()
    return low.startswith("ami_") or ("ami" in low and "_to_" not in low)


def _is_mapped_wall_patch(name: str, patch_types: dict[str, str] | None) -> bool:
    if patch_types and patch_types.get(name) == "mappedWall":
        return True
    return "_to_" in name


# Geometric constraint patch types: the field entry MUST repeat the patch
# type or OpenFOAM aborts at startup with a type mismatch.
_CONSTRAINT_PATCH_TYPES = ("symmetryPlane", "empty", "wedge", "cyclic")


def _constraint_rule(name: str, patch_types: dict[str, str] | None) -> str | None:
    """Field rule for constraint patches (``symmetryPlane``/``empty``/...)."""
    if not patch_types:
        return None
    ptype = patch_types.get(name)
    if ptype in _CONSTRAINT_PATCH_TYPES:
        return f"type {ptype};"
    return None


def _is_opening_patch(name: str, opening_patches: set[str] | None) -> bool:
    """Total-pressure opening: ``open*`` name or CGNS inlet/outlet BC."""
    if opening_patches and name in opening_patches:
        return True
    return name == "open" or name.startswith("open")


def _field_T(
    region_type: str,
    patches: list[str],
    T0: float = 300.0,
    *,
    patch_types: dict[str, str] | None = None,
    opening_patches: set[str] | None = None,
    convection_patches: dict[str, tuple[float, float]] | None = None,
) -> str:
    """0/<region>/T.

    *convection_patches* maps patch name → ``(Ta, h)`` and selects
    ``externalWallHeatFluxTemperature`` (mode ``coefficient``) for outer
    walls cooled by ambient convection.
    """
    kappa = "fluidThermo" if region_type == "fluid" else "solidThermo"
    convection_patches = convection_patches or {}
    blocks: list[str] = []
    for p in patches:
        rule = _constraint_rule(p, patch_types)
        if rule is not None:
            blocks.append(f"    {p}\n    {{\n        {rule}\n    }}")
        elif _is_cyclic_ami_patch(p, patch_types):
            blocks.append(f"    {p}\n    {{\n        type            cyclicAMI;\n    }}")
        elif _is_mapped_wall_patch(p, patch_types):
            blocks.append(
                f"""    {p}
    {{
        type            compressible::turbulentTemperatureRadCoupledMixed;
        Tnbr            T;
        kappaMethod     {kappa};
        value           $internalField;
    }}"""
            )
        elif p in convection_patches:
            ta, htc = convection_patches[p]
            blocks.append(
                f"""    {p}
    {{
        type            externalWallHeatFluxTemperature;
        mode            coefficient;
        Ta              constant {ta};
        h               uniform {htc};
        kappaMethod     {kappa};
        qr              none;
        value           $internalField;
    }}"""
            )
        elif _is_opening_patch(p, opening_patches):
            blocks.append(
                f"""    {p}
    {{
        type            inletOutlet;
        inletValue      uniform {T0};
        value           uniform {T0};
    }}"""
            )
        else:
            blocks.append(
                f"""    {p}
    {{
        type            zeroGradient;
    }}"""
            )
    body = '    #includeEtc "caseDicts/setConstraintTypes"\n\n' + (
        "\n\n".join(blocks) if blocks else "    // (no patches yet)"
    )
    return (
        _foam_header("volScalarField", "T", "0")
        + f"""
dimensions      [0 0 0 1 0 0 0];
internalField   uniform {T0};

boundaryField
{{
{body}
}}

// ************************************************************************* //
"""
    )


def _field_U(
    patches: list[str],
    *,
    patch_types: dict[str, str] | None = None,
    moving_wall_patches: list[str] | None = None,
    opening_patches: set[str] | None = None,
) -> str:
    moving = set(moving_wall_patches or [])
    blocks: list[str] = []
    for p in patches:
        rule = _constraint_rule(p, patch_types)
        if rule is not None:
            blocks.append(f"    {p}\n    {{\n        {rule}\n    }}")
        elif _is_cyclic_ami_patch(p, patch_types):
            blocks.append(f"    {p}\n    {{\n        type            cyclicAMI;\n    }}")
        elif _is_mapped_wall_patch(p, patch_types):
            blocks.append(
                f"""    {p}
    {{
        type            noSlip;
    }}"""
            )
        elif _is_opening_patch(p, opening_patches):
            blocks.append(
                f"""    {p}
    {{
        type            pressureInletOutletVelocity;
        value           uniform (0 0 0);
    }}"""
            )
        elif p in moving or "impeller" in p.lower():
            # MRF: blade walls use absolute movingWallVelocity (omega x r)
            blocks.append(
                f"""    {p}
    {{
        type            movingWallVelocity;
        value           uniform (0 0 0);
    }}"""
            )
        else:
            blocks.append(
                f"""    {p}
    {{
        type            noSlip;
    }}"""
            )
    body = '    #includeEtc "caseDicts/setConstraintTypes"\n\n' + (
        "\n\n".join(blocks) if blocks else "    // (no patches yet)"
    )
    return (
        _foam_header("volVectorField", "U", "0")
        + f"""
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);

boundaryField
{{
{body}
}}

// ************************************************************************* //
"""
    )


def _field_p(
    patches: list[str],
    p0: float = 101325.0,
    *,
    patch_types: dict[str, str] | None = None,
    opening_patches: set[str] | None = None,
) -> str:
    blocks: list[str] = []
    for p in patches:
        rule = _constraint_rule(p, patch_types)
        if rule is not None:
            blocks.append(f"    {p}\n    {{\n        {rule}\n    }}")
        elif _is_cyclic_ami_patch(p, patch_types):
            blocks.append(f"    {p}\n    {{\n        type            cyclicAMI;\n    }}")
        elif _is_opening_patch(p, opening_patches):
            blocks.append(
                f"""    {p}
    {{
        type            calculated;
        value           uniform {p0};
    }}"""
            )
        else:
            blocks.append(
                f"""    {p}
    {{
        type            calculated;
        value           $internalField;
    }}"""
            )
    body = '    #includeEtc "caseDicts/setConstraintTypes"\n\n' + (
        "\n\n".join(blocks) if blocks else "    // (no patches)"
    )
    return (
        _foam_header("volScalarField", "p", "0")
        + f"""
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform {p0};

boundaryField
{{
{body}
}}

// ************************************************************************* //
"""
    )


def _field_p_rgh(
    patches: list[str],
    p0: float = 101325.0,
    *,
    patch_types: dict[str, str] | None = None,
    opening_patches: set[str] | None = None,
) -> str:
    blocks: list[str] = []
    for p in patches:
        rule = _constraint_rule(p, patch_types)
        if rule is not None:
            blocks.append(f"    {p}\n    {{\n        {rule}\n    }}")
        elif _is_cyclic_ami_patch(p, patch_types):
            blocks.append(f"    {p}\n    {{\n        type            cyclicAMI;\n    }}")
        elif _is_opening_patch(p, opening_patches):
            blocks.append(
                f"""    {p}
    {{
        type            prghTotalPressure;
        p0              uniform {p0};
        value           uniform {p0};
    }}"""
            )
        else:
            blocks.append(
                f"""    {p}
    {{
        type            fixedFluxPressure;
        value           $internalField;
    }}"""
            )
    body = '    #includeEtc "caseDicts/setConstraintTypes"\n\n' + (
        "\n\n".join(blocks) if blocks else "    // (no patches)"
    )
    return (
        _foam_header("volScalarField", "p_rgh", "0")
        + f"""
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform {p0};

boundaryField
{{
{body}
}}

// ************************************************************************* //
"""
    )


# ---------------------------------------------------------------------------
# Allrun scripts
# ---------------------------------------------------------------------------


def _decompose_par_dict(n_procs: int = 8, *, location: str = "system") -> str:
    return (
        _foam_header("dictionary", "decomposeParDict", location)
        + f"""
numberOfSubdomains  {n_procs};

method          scotch;

// ************************************************************************* //
"""
    )


def _allclean() -> str:
    return """#!/bin/sh
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions
cleanCase
#------------------------------------------------------------------------------
"""


# ---------------------------------------------------------------------------
# Public API
