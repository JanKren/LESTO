# LESTO – OpenFOAM / GEMS Coupling Development

This repository contains the OpenFOAM cases and development steps used for
the LESTO project.

The long-term goal is to couple OpenFOAM with GEMS/GEMS3K for simulations
involving fluid flow together with thermodynamic equilibrium calculations.

The current pipe-flow cases are deliberately simple. They are evolutionary
steps used to establish and validate the OpenFOAM side of the workflow before
introducing the GEMS coupling and the full LESTO physics.

## Development strategy

The repository is organized as a sequence of increasingly complex cases.

Current steps:

- `laminarPipeFlow001`  
  Basic steady laminar pipe flow with uniform inlet and uniform initial field.

- `laminarPipeFlow002`  
  Same flow problem with a parabolic inlet velocity profile.

- `laminarPipeFlow003`  
  Parabolic inlet and parabolic initial velocity field.

Further cases will progressively introduce the additional physics and finally
the OpenFOAM–GEMS coupling required for the LESTO project.

## Repository structure

    run/
        laminarPipeFlow001/
        laminarPipeFlow002/
        laminarPipeFlow003/
        ...

Each case contains its own `readme.md` describing what changed relative to the
previous development step.

## Typical OpenFOAM workflow

From inside a case directory:

    gmsh -3 pipe.geo
    gmshToFoam pipe.msh
    checkMesh
    simpleFoam > out &

To follow the run:

    tail -f out

## Generated files

Generated OpenFOAM data are not tracked, including:

    constant/polyMesh/
    numerical time directories
    dynamicCode/
    *.foam

These can be regenerated from the committed case input files.

Solver output files such as `out` are kept as a record of the actual run.

## Software

Current development uses:

    OpenFOAM v2606
    Gmsh
    ParaView

Later stages will also include:

    GEMS / GEMS3K

The phase-change solver's M10a stage now supports independent PbI2/BiI3
pairs, shared multi-species sources and element conservation outputs.
See [the solver documentation](applications/rhoFixedFlowFoam/readme.md)
and [the preliminary channel](run/030-pbi2-bii3-channel/readme.md).
M10b adds bridge-free homogeneous gas re-speciation and reaction/element
accounts; see [run/031](run/031-pbbi-gas-speciation/readme.md).
M10c adds optional frozen/local GEMS gas engines; see
[run/032](run/032-pbbi-gems-speciation/readme.md). Reactive/shared-condensate
coupling remains for M10d before experimental validation.
