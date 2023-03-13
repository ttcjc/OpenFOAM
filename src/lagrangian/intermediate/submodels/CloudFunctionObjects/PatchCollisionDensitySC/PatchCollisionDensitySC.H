/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::PatchCollisionDensitySC

Description
    Function object which generates fields of the number and rate of collisions
    per unit area on all patches. Can optionally take a minimum speed below
    which a collision is not counted.

    Example usage:
    \verbatim
    PatchCollisionDensitySC1
    {
        type        PatchCollisionDensitySC;
        minSpeed    1e-3;
    }
    \endverbatim

SourceFiles
    PatchCollisionDensitySC.C

\*---------------------------------------------------------------------------*/

#ifndef PatchCollisionDensitySC_H
#define PatchCollisionDensitySC_H

#include "CloudFunctionObjectSC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PatchCollisionDensitySC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PatchCollisionDensitySC
:
    public CloudFunctionObjectSC<CloudType>
{
    // Private Data

        typedef typename CloudType::particleType parcelType;

        //- The threshold for a collision
        const scalar minSpeed_;

        //- The field of the number of collisions per unit area
        volScalarField::Boundary collisionDensity_;

        //- The field of the number of collisions per unit area at the last
        //  output
        volScalarField::Boundary collisionDensity0_;

        //- The time at the last output
        scalar time0_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("PatchCollisionDensitySC");


    // Constructors

        //- Construct from dictionary
        PatchCollisionDensitySC
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        PatchCollisionDensitySC(const PatchCollisionDensitySC<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObjectSC<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObjectSC<CloudType>>
            (
                new PatchCollisionDensitySC<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PatchCollisionDensitySC();


    // Member Functions

        // Evaluation

            //- Post-patch hook
            virtual void postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchCollisionDensitySC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
