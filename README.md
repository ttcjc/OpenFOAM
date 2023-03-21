# OpenFOAM

Modifications to OpenFOAM v7, improving Lagrangian multi-phase functionality for use in the study of automotive surface contamination at Loughborough University.

#### Existing OpenFOAM Files

Any minor additions or alterations made to existing OpenFOAM files are denoted by the inclusion of *CJC* comments.

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

Brings DESModelRegions functionality from OpenFOAM v2.3 to OpenFOAM v7 based on [wyldckat's](https://github.com/wyldckat/DESModelRegions) own implementation for OpenFOAM v5 and OpenFOAM v6.

Outputs an indicator of RAS and LES regions within a DES simulation:
- '0' indicates RAS
- '1' indicates LES

##### forceCoeffsExtended

Expands forceCoeff functionality to include outputs for all three force coefficients (lift, drag and side force) as well as all three moment coefficients (yaw, pitch and roll)

### lagrangian

#### 'SC' (Surface Contamination) File Variants

Variations of existing OpenFOAM files developed to support custom solver (pisoFoamSC) functionality without impacting native operation.

Additional functionality includes:
- Ability to record Lagrangian data at user-defined time interval (separate from Eulerian phase)
- Ability to record Lagrangian data to separate file when particles impinge on a surface in a user-defined volume (and subsequent removal of particles from simulation)
- Ability to record Lagrangian data to separate file(s) as particles cross user-defined streamwise plane(s) of interest (and subsequent removal of particles from simulation)

##### WheelInjectionSC

A custom injector that introduces particles randomly about a user-defined cylinder representing a wheel.

Two injection methods are available, intended to produce injection patterns representative of both tread pickup and capillary action, as defined by Weir et al. in "Reduction of Adverse Aerodynamic Effects of Large Trucks"
