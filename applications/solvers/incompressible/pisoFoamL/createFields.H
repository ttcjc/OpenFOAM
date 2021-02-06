/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

// Pressure
Info << "Reading Field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Velocity
Info << "Reading Field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Read transportProperties Dict
Info << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
		"transportProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

// Store Fluid Density
dimensionedScalar rhoValue
(
	transportProperties.lookup("rhoInf")
);

// Store Reciprocal of Density (Two-Way Coupling)
dimensionedScalar rhoInverse
(
	"rhoInverse",
	(1/rhoValue)
);

// Create Field for Fluid Density
volScalarField rho
(
    IOobject
    (
		"rho",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	rhoValue
);

// Create Field for Dynamic Viscosity
volScalarField mu
(
    IOobject
    (
		"mu",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	laminarTransport.nu()*rhoValue
);

#include "readGravitationalAcceleration.H"

// Initialise Particle Cloud
const word kinematicCloudName
(
	args.optionLookupOrDefault<word>("cloudName", "kinematicCloud")
);

Info << "Initialising Particle Cloud " << kinematicCloudName << nl << endl;

basicKinematicCollidingCloud kinematicCloud
(
	kinematicCloudName,
	rho,
	U,
	mu,
	g
);

// ************************************************************************* //
