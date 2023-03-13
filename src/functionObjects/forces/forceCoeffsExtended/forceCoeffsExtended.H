/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::forceCoeffsExtended

Description
    Extends the forces functionObject by providing lift, drag and moment
    coefficients.  The data can optionally be output into bins, defined in a
    given direction.

    Example of function object specification:
    \verbatim
    forceCoeffsExtended
    {
        type        forceCoeffsExtended;
        libs        ("libforcesExtended.so");
        ...
        log         yes;
        patches     (walls);
        liftDir     (0 0 1);
        dragDir     (1 0 0);
        sideDir_    (0 1 0);
        pitchAxis   (0 1 0);
        yawAxis     (0 0 1);
        rollAxis    (1 0 0);
        magUInf     1;
        lRef        1;
        Aref        1;
    }
    \endverbatim

Usage
    \table
        Property    | Description                                   | Required  | Default value
        type        | type name: forceCoeffsExtended                | yes       |
        log         | write force data to standard output           | no        | no
        patches     | patches included in the forces calculation    | yes       |
        liftDir     | lift direction                                | yes       | (0 0 1)
        dragDir     | drag direction                                | yes       | (1 0 0)
        sideDir     | side force direction                          | yes       | (0 1 0)
        pitchAxis   | pitch axis                                    | yes       | (0 1 0)
        yawAxis     | yaw axis                                      | yes       | (0 0 1)
        rollAxis    | roll axis                                     | yes       | (1 0 0)
        magUInf     | freestream velocity magnitude                 | yes       | 1
        lRef        | reference length for moment calculations      | yes       | 1
        Aref        | reference area                                | yes       | 1
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::timeControl
    Foam::functionObjects::forces

SourceFiles
    forceCoeffsExtended.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forceCoeffsExtended_H
#define functionObjects_forceCoeffsExtended_H

#include "forces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class forceCoeffsExtended Declaration
\*---------------------------------------------------------------------------*/

class forceCoeffsExtended
:
    public forces
{
    // Private Data

        // Force coefficient geometry

            //- Lift
            vector liftDir_;

            //- Drag
            vector  dragDir_;

            //- Side Force
            vector  sideDir_;

            //- Pitch
            vector pitchAxis_;

            //- Yaw
            vector yawAxis_;

            //- Roll
            vector rollAxis_;


        // Freestream conditions

            //- Velocity magnitude
            scalar magUInf_;


        // Reference scales

            //- Length
            scalar lRef_;

            //- Area
            scalar Aref_;


protected:

    //- Output file header information
    virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("forceCoeffsExtended");


    // Constructors

        //- Construct from Time and dictionary
        forceCoeffsExtended
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        forceCoeffsExtended(const forceCoeffsExtended&) = delete;


    //- Destructor
    virtual ~forceCoeffsExtended();


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the forces
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const forceCoeffsExtended&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
