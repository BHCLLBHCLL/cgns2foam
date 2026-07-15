"""Generate chtMultiRegionSimpleFoam case scaffolding from a mesh + coupling scan.

Writes regionProperties, per-region thermo/system/0.orig stubs, CHT controlDict,
AMI createPatchDict (when needed), and Allrun / Allrun.pre / Allclean scripts.
Mesh prep (``splitMeshRegions`` / optional ``createPatch``) is deferred to
Allrun.pre inside an OpenFOAM v2412 environment.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .couplings import CouplingMethod, CouplingReport
from .topology import Mesh


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


def _thermophysical_fluid() -> str:
    t = _DEFAULT_FLUID_MAT["thermoType"]
    m = _DEFAULT_FLUID_MAT["mixture"]
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


def _thermophysical_solid() -> str:
    t = _DEFAULT_SOLID_MAT["thermoType"]
    m = _DEFAULT_SOLID_MAT["mixture"]
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


def _control_dict_cht(*, end_time: int = 500, write_interval: int = 50) -> str:
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
purgeWrite      0;

writeFormat     binary;
writePrecision  8;
writeCompression off;

timeFormat      general;
timePrecision   8;
runTimeModifiable true;

functions
{{
}}

// ************************************************************************* //
"""
    )


def _fv_schemes_fluid() -> str:
    return (
        _foam_header("dictionary", "fvSchemes", "system")
        + """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
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
laplacianSchemes { default Gauss linear limited 0.33; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.33; }

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
    nNonOrthogonalCorrectors 0;
    residualControl { default 1e-5; }
}
relaxationFactors
{
    fields { p_rgh 0.3; }
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
    residualControl { default 1e-5; }
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


def _field_T(region_type: str, patches: list[str], T0: float = 300.0) -> str:
    kappa = "fluidThermo" if region_type == "fluid" else "solidThermo"
    blocks: list[str] = []
    for p in patches:
        if "ami" in p.lower() and "_to_" not in p:
            blocks.append(f"    {p}\n    {{\n        type            cyclicAMI;\n    }}")
        elif "_to_" in p:
            blocks.append(
                f"""    {p}
    {{
        type            compressible::turbulentTemperatureRadCoupledMixed;
        Tnbr            T;
        kappaMethod     {kappa};
        value           $internalField;
    }}"""
            )
        else:
            blocks.append(
                f"""    {p}
    {{
        type            zeroGradient;
    }}"""
            )
    body = "\n\n".join(blocks) if blocks else "    // (no patches yet)"
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


def _field_U(patches: list[str]) -> str:
    blocks: list[str] = []
    for p in patches:
        if "ami" in p.lower() and "_to_" not in p:
            blocks.append(f"    {p}\n    {{\n        type            cyclicAMI;\n    }}")
        elif "_to_" in p:
            blocks.append(
                f"""    {p}
    {{
        type            noSlip;
    }}"""
            )
        elif p.startswith("open"):
            blocks.append(
                f"""    {p}
    {{
        type            pressureInletOutletVelocity;
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
    body = "\n\n".join(blocks) if blocks else "    // (no patches yet)"
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


def _field_p(patches: list[str], p0: float = 101325.0) -> str:
    blocks: list[str] = []
    for p in patches:
        if "ami" in p.lower() and "_to_" not in p:
            blocks.append(f"    {p}\n    {{\n        type            cyclicAMI;\n    }}")
        else:
            blocks.append(
                f"""    {p}
    {{
        type            calculated;
        value           $internalField;
    }}"""
            )
    body = "\n\n".join(blocks) if blocks else "    // (no patches)"
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


def _field_p_rgh(patches: list[str]) -> str:
    blocks: list[str] = []
    for p in patches:
        if "ami" in p.lower() and "_to_" not in p:
            blocks.append(f"    {p}\n    {{\n        type            cyclicAMI;\n    }}")
        elif p.startswith("open"):
            blocks.append(
                f"""    {p}
    {{
        type            fixedValue;
        value           uniform 0;
    }}"""
            )
        else:
            blocks.append(
                f"""    {p}
    {{
        type            fixedFluxPressure;
        value           uniform 0;
    }}"""
            )
    body = "\n\n".join(blocks) if blocks else "    // (no patches)"
    return (
        _foam_header("volScalarField", "p_rgh", "0")
        + f"""
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 0;

boundaryField
{{
{body}
}}

// ************************************************************************* //
"""
    )


def _create_patch_ami(
    pairs: list[tuple[str, str]],
    axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
) -> str:
    blocks: list[str] = []
    for master, slave in pairs:
        for name, nbr in ((master, slave), (slave, master)):
            blocks.append(
                f"""    {{
        name {name};
        patchInfo
        {{
            type            cyclicAMI;
            neighbourPatch  {nbr};
            transform       noOrdering;
            matchTolerance  0.001;
            rotationAxis    ({axis[0]} {axis[1]} {axis[2]});
        }}
        constructFrom patches;
        patches ({name});
    }}"""
            )
    body = "\n".join(blocks)
    return (
        _foam_header("dictionary", "createPatchDict")
        + f"""
pointSync false;

patches
(
{body}
);

// ************************************************************************* //
"""
    )


# ---------------------------------------------------------------------------
# Allrun scripts
# ---------------------------------------------------------------------------


def _allrun_pre_full(
    fluid: list[str],
    solid: list[str],
    *,
    has_ami: bool,
) -> str:
    region_names = fluid + solid
    lines = [
        "#!/bin/sh",
        "set -e",
        'cd "${0%/*}" || exit',
        ". ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions",
        "#------------------------------------------------------------------------------",
        "# cgns2foam – CHT mesh prep for chtMultiRegionSimpleFoam",
        "",
    ]
    if has_ami:
        lines.append("runApplication -s createPatch createPatch -overwrite")
    lines.extend(
        [
            "cp -f system/regionProperties constant/regionProperties",
            "runApplication -s splitMeshRegions splitMeshRegions -cellZonesOnly -overwrite",
            "",
            "# Deploy per-region constant/system and CHT controlDict",
            f"for region in {' '.join(region_names)}; do",
            '    mkdir -p "constant/${region}" "system/${region}"',
            '    cp -f constant.orig/"${region}"/* "constant/${region}/" 2>/dev/null || true',
            '    cp -f system.orig/"${region}"/* "system/${region}/" 2>/dev/null || true',
            "done",
            "cp -f system/controlDict.cht system/controlDict",
            "",
            "restore0Dir -allRegions",
            "",
        ]
    )
    for reg in solid:
        for f in ("U", "p_rgh", "k", "epsilon", "nut", "alphat"):
            lines.append(f'rm -f "0/{reg}/{f}" 2>/dev/null || true')
    lines.append("")
    lines.append("#------------------------------------------------------------------------------")
    return "\n".join(lines) + "\n"


def _allrun() -> str:
    return """#!/bin/sh
set -e
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions
#------------------------------------------------------------------------------
./Allrun.pre
runApplication $(getApplication)
#------------------------------------------------------------------------------
"""


def _allclean() -> str:
    return """#!/bin/sh
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions
cleanCase
#------------------------------------------------------------------------------
"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def write_cht_case(
    out_dir: str | Path,
    mesh: Mesh,
    report: CouplingReport,
    *,
    end_time: int = 500,
    write_interval: int = 50,
    gravity: list[float] | None = None,
) -> dict[str, Any]:
    """Augment an already-written mono-block OpenFOAM case with CHT scaffolding.

    Expects ``constant/polyMesh`` already present (from :func:`write_case`).
    """
    out = Path(out_dir)
    fluid = list(report.fluid_regions)
    solid = list(report.solid_regions)
    if not fluid and not solid:
        raise ValueError("Coupling report has no fluid/solid regions")

    region_by_foam = {r.foam_name: r for r in report.regions}
    patch_names = [p.name for p in mesh.patches]

    # Map foam region → patches that belong to that zone's original BCs
    # (approximate: patches whose name stem matches a BC of that region)
    def patches_for_region(foam_name: str) -> list[str]:
        reg = region_by_foam.get(foam_name)
        if not reg:
            return list(patch_names)
        bc_set = set(reg.bc_names)
        owned = [p for p in patch_names if p in bc_set or any(
            p == b or p.startswith(b + "_") for b in bc_set
        )]
        return owned if owned else list(patch_names)

    # --- top-level system / constant ---
    _write_text(out / "system" / "regionProperties", _region_properties(fluid, solid))
    _write_text(out / "system" / "controlDict.cht", _control_dict_cht(
        end_time=end_time, write_interval=write_interval,
    ))
    # Keep monolithic controlDict usable by mesh utilities until Allrun.pre
    # replaces it with controlDict.cht after split.
    _write_text(out / "constant" / "g", _gravity(gravity))

    ami_pairs: list[tuple[str, str]] = []
    for c in report.couplings:
        if c.method != CouplingMethod.CYCLIC_AMI:
            continue
        # Patch names in mono mesh are disambiguated; find BC-based names
        master_patches = [
            p.name for p in mesh.patches
            if p.name == c.master_bc or p.name.startswith(c.master_bc + "_")
        ]
        slave_patches = [
            p.name for p in mesh.patches
            if p.name == c.slave_bc or p.name.startswith(c.slave_bc + "_")
        ]
        if master_patches and slave_patches:
            ami_pairs.append((master_patches[0], slave_patches[0]))
        else:
            # Fall back to raw BC names (createPatch may still find them)
            ami_pairs.append((c.master_bc, c.slave_bc))

    # Dedupe AMI pairs
    seen_ami: set[frozenset[str]] = set()
    unique_ami: list[tuple[str, str]] = []
    for a, b in ami_pairs:
        key = frozenset({a, b})
        if key in seen_ami:
            continue
        seen_ami.add(key)
        unique_ami.append((a, b))

    if unique_ami:
        _write_text(out / "system" / "createPatchDict", _create_patch_ami(unique_ami))

    # --- per-region staged files ---
    for foam_name in fluid + solid:
        reg = region_by_foam[foam_name]
        rtype = reg.region_type
        cdir = out / "constant.orig" / foam_name
        sdir = out / "system.orig" / foam_name
        odir = out / "0.orig" / foam_name
        owned = patches_for_region(foam_name)

        if rtype == "fluid":
            _write_text(cdir / "thermophysicalProperties", _thermophysical_fluid())
            _write_text(cdir / "turbulenceProperties", _turbulence_laminar())
            _write_text(sdir / "fvSchemes", _fv_schemes_fluid())
            _write_text(sdir / "fvSolution", _fv_solution_fluid())
            _write_text(odir / "T", _field_T("fluid", owned))
            _write_text(odir / "U", _field_U(owned))
            _write_text(odir / "p", _field_p(owned))
            _write_text(odir / "p_rgh", _field_p_rgh(owned))
        else:
            _write_text(cdir / "thermophysicalProperties", _thermophysical_solid())
            _write_text(sdir / "fvSchemes", _fv_schemes_solid())
            _write_text(sdir / "fvSolution", _fv_solution_solid())
            _write_text(odir / "T", _field_T("solid", owned))
            _write_text(odir / "p", _field_p(owned))

    summary = report.to_dict()
    summary["ami_pairs"] = [{"master": a, "slave": b} for a, b in unique_ami]
    summary["solver"] = "chtMultiRegionSimpleFoam"
    _write_text(out / "setup_report.json", json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    _write_text(out / "coupling_scan.json", json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    _write_text(
        out / "Allrun.pre",
        _allrun_pre_full(fluid, solid, has_ami=bool(unique_ami)),
    )
    _write_text(out / "Allrun", _allrun())
    _write_text(out / "Allclean", _allclean())

    # Make scripts executable on POSIX (no-op on Windows)
    for name in ("Allrun", "Allrun.pre", "Allclean"):
        path = out / name
        try:
            mode = path.stat().st_mode
            path.chmod(mode | 0o111)
        except OSError:
            pass

    return summary
