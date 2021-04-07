/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "LESresolution.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(LESresolution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        LESresolution,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::LESresolution::calc()
{
    if(foundObject<volSymmTensorField>(fieldName_))
    {
        if
        (
            foundObject<volSymmTensorField>("UPrime2Mean")
            &&
            foundObject<volSymmTensorField>("R")
        )
        {
            const volSymmTensorField& UPrime2Mean = lookupObject<volSymmTensorField>("UPrime2Mean");
            const volSymmTensorField& R = lookupObject<volSymmTensorField>("R");

            return store
            (
                resultName_,
                tr(UPrime2Mean) / max(tr(UPrime2Mean + R), dimensionedScalar(dimensionSet(0, 2, -2, 0, 0, 0, 0), vSmall))
            );
        }

        else
        {
            return false;
        }
    }

    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::LESresolution::LESresolution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "LESresolution")
{
    setResultName(typeName, "LESresolution");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::LESresolution::~LESresolution()
{}


// ************************************************************************* //
