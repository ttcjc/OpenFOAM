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

\*---------------------------------------------------------------------------*/

#include "RelativeVelocitySC.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::RelativeVelocitySC<CloudType>::write()
{
    IOField<vector> URel
    (
        this->owner().fieldIOobject("URel", IOobject::NO_READ),
        this->owner().size()
    );

    autoPtr<interpolation<vector>> UcInterp
    (
        interpolation<vector>::New
        (
            this->owner().solution().interpolationSchemes(),
            this->owner().U()
        )
    );

    label i = 0;
    forAllConstIter(typename CloudType, this->owner(), iter)
    {
        URel[i] =
            iter().U()
          - UcInterp->interpolate
            (
                iter().coordinates(),
                iter().currentTetIndices()
            );
        ++ i;
    }

    URel.write(this->owner().size() > 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RelativeVelocitySC<CloudType>::RelativeVelocitySC
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObjectSC<CloudType>(dict, owner, modelName, typeName)
{}


template<class CloudType>
Foam::RelativeVelocitySC<CloudType>::RelativeVelocitySC
(
    const RelativeVelocitySC<CloudType>& ppm
)
:
    CloudFunctionObjectSC<CloudType>(ppm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RelativeVelocitySC<CloudType>::~RelativeVelocitySC()
{}


// ************************************************************************* //
