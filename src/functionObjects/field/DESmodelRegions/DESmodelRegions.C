/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "DESmodelRegions.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentFluidThermoModels.H"
#include "wallDist.H"
#include "bound.H"
#include "fvcGrad.H"
#include "SpalartAllmarasDES.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(DESmodelRegions, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        DESmodelRegions,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::DESmodelRegions::DESmodelRegions
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::DESmodelRegions::~DESmodelRegions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::functionObjects::DESmodelRegions::modelName()
{
    return Foam::turbulenceModel::propertiesName;
}


bool Foam::functionObjects::DESmodelRegions::execute()
{
    if
    (
        obr_.foundObject<incompressible::turbulenceModel>(modelName())
    )
    {
        const incompressible::turbulenceModel& turbModel =
        obr_.lookupObject<incompressible::turbulenceModel>(modelName());

        typedef LESModels::SpalartAllmarasDES<incompressible::turbulenceModel> icoDESmodel;

        if
        (
            isA<icoDESmodel>(turbModel)
        )
        {
            const icoDESmodel& DES = dynamic_cast<const icoDESmodel&>(turbModel);

            tmp<volScalarField> tDESfield(DES.LESRegion());
            volScalarField &DESmodelRegions = tDESfield.ref();

            word result("DESmodelRegions");

            return store
            (
                result,
                DESmodelRegions / dimensionedScalar(dimless, 1)
            );
        }
    }

    else if
    (
        obr_.foundObject<compressible::turbulenceModel>(modelName())
    )
    {
        FatalErrorInFunction
            << "Support for Compressible Solvers Not Yet Implemented"
            << exit(FatalError);
    }

    else
    {
        FatalErrorInFunction
            << "Invalid Turbulence Model"
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::DESmodelRegions::write()
{
    Log << type() << " " << name() << " Output:" << nl
        << "    Writing DES Field" << nl << endl;

    const volScalarField& DESmodelRegions = obr_.lookupObject<volScalarField>("DESmodelRegions");
    DESmodelRegions.write();

    return true;
}

// ************************************************************************* //
