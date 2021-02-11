/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
	Foam::functionObjects::LESresolution

Description
	Outputs the ratio of resolved turbulent kinetic energy to total
	turbulent kinetic energy within LES or DES simulations

	Dependencies:
	- fieldAverage function to obtain UPrime2Mean (Resolved Reynolds Stress Tensor)
	- turbulenceFields function to obtain R (Subgrid Reynolds Stress Tensor)

See also
	Foam::functionObjects::fieldExpression
	Foam::functionObjects::fvMeshFunctionObject

SourceFiles
	LESresolution.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_LESresolution_H
#define functionObjects_LESresolution_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class LESresolution Declaration
\*---------------------------------------------------------------------------*/

class LESresolution
:
	public fieldExpression
{
	// Private Member Functions

		//- Calculate the LESresolution field and return true if successful
		virtual bool calc();


public:

	//- Runtime type information
	TypeName("LESresolution");


	// Constructors

		//- Construct from Time and dictionary
		LESresolution
		(
			const word& name,
			const Time& runTime,
			const dictionary& dict
		);


	// Destructor
	virtual ~LESresolution();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
