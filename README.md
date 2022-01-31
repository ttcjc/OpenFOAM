# OpenFOAM

Modifications to OpenFOAM v7, improving Lagrangian multi-phase functionality for use in the study of automotive surface contamination at Loughborough University

#### Existing OpenFOAM Files

Any minor additions or alterations made to existing OpenFOAM files are denoted by the inclusion of *CJC* comments

For example:

```
This is a single line modification // CJC
```

```
// CJC {
    This is a multi-line modification #1
    This is a multi-line modification #2
// } CJC
```

## applications

### solvers

##### pisoFoamSC

Transient solver for incompressible, particle-laden, turbulent flow,
using the PISO algorithm.

Sub-models include:
- Turbulence modelling, i.e. laminar, RAS or LES
- Run-time selectable MRF and finite volume options, e.g. explicit porosity
- Optional modelling of Lagrangian particle clouds

## src

### functionObjects

##### DESModelRegions

Brings DESModelRegions functionality from OpenFOAM v2.3 to OpenFOAM v7 based on [wyldckat's](https://github.com/wyldckat/DESModelRegions) own implementation for OpenFOAM v5 and OpenFOAM v6

Outputs an indicator of RAS and LES regions within a DES simulation:
- '0' indicates RAS
- '1' indicates LES

##### LESresolution

Outputs the ratio of resolved turbulent kinetic energy to total turbulent kinetic energy within LES or DES simulations

Dependencies:
- fieldAverage function to obtain UPrime2Mean (Resolved Reynolds Stress Tensor)
- turbulenceFields function to obtain R (Subgrid Reynolds Stress Tensor)

### lagrangian

#### 'SC' (Surface Contamination) File Variants

Variations of existing OpenFOAM files developed to support custom solver (pisoFoamSC) functionality without impacting native operation
