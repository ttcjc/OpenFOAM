/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::functionObjects::DESmodelRegions

Description
    Brings DESModelRegions functionality from OpenFOAM v2.3 to OpenFOAM v7
    based on wyldckat's own implementation for OpenFOAM v5 and OpenFOAM v6

    Outputs an indicator of RAS and LES regions within a DES simulation:
    - '0' indicates RAS
    - '1' indicates LES

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    DESmodelRegions.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_DESmodelRegions_H
#define functionObjects_DESmodelRegions_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class DESmodelRegions Declaration
\*---------------------------------------------------------------------------*/

class DESmodelRegions
:
    public fvMeshFunctionObject
{
public:

    static const word& modelName();

    //- Runtime type information
    TypeName("DESmodelRegions");


    // Constructors

        //- Construct from Time and dictionary
        DESmodelRegions
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~DESmodelRegions();


    // Member Functions

        //- Calculate DES field
        virtual bool execute();

        //- Write DES field
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator = (const DESmodelRegions&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
