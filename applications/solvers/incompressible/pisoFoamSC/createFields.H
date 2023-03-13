/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

Info<< "Reading field p\n" << endl;
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

Info<< "Reading field U\n" << endl;
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

// CJC {
    // Read transportProperties
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

    // Store value of carrier fluid density from transportProperties
    dimensionedScalar rhoValue
    (
        transportProperties.lookup("rhoInf")
    );

    // Store reciprocal of carrier fluid density (used for two-wap coupling)
    dimensionedScalar rhoInverseValue
    (
        "rhoInverseValue",
        (1 / rhoValue)
    );

    // Create volScalarField of carrier fluid density
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

    // Create volScalarField for carrier fluid dynamic viscosity
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

    // Allow Lagrangian particles to be influenced by gravity
    #include "readGravitationalAcceleration.H"

    // Initialise Lagrangian particle cloud
    const word kinematicCloudName
    (
        args.optionLookupOrDefault<word>("cloudName", "kinematicCloud")
    );

    Info << "Initialising Lagrangian Particle Cloud " << kinematicCloudName << nl << endl;

    basicKinematicCloudSC kinematicCloud
    (
        kinematicCloudName,
        rho,
        U,
        mu,
        g
    );
// } CJC

// ************************************************************************* //
