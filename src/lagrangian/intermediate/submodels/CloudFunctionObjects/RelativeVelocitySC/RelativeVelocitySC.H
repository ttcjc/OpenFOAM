/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::RelativeVelocitySC

Description
    Function object which generates a lagrangian field containing the velocity
    of each particle relative to the velocity of the surrounding fluid.

    Example usage:
    \verbatim
    RelativeVelocitySC1
    {
        type        RelativeVelocitySC;
    }
    \endverbatim

SourceFiles
    RelativeVelocitySC.C

\*---------------------------------------------------------------------------*/

#ifndef RelativeVelocitySC_H
#define RelativeVelocitySC_H

#include "CloudFunctionObjectSC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class RelativeVelocitySC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RelativeVelocitySC
:
    public CloudFunctionObjectSC<CloudType>
{
protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("RelativeVelocitySC");


    // Constructors

        //- Construct from dictionary
        RelativeVelocitySC
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        RelativeVelocitySC(const RelativeVelocitySC<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObjectSC<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObjectSC<CloudType>>
            (
                new RelativeVelocitySC<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RelativeVelocitySC();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RelativeVelocitySC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
