/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "VoidFractionSC.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFractionSC<CloudType>::write()
{
    if (thetaPtr_.valid())
    {
        thetaPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "thetaPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFractionSC<CloudType>::VoidFractionSC
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObjectSC<CloudType>(dict, owner, modelName, typeName),
    thetaPtr_(nullptr)
{}


template<class CloudType>
Foam::VoidFractionSC<CloudType>::VoidFractionSC
(
    const VoidFractionSC<CloudType>& vf
)
:
    CloudFunctionObjectSC<CloudType>(vf),
    thetaPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFractionSC<CloudType>::~VoidFractionSC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFractionSC<CloudType>::preEvolve()
{
    if (thetaPtr_.valid())
    {
        thetaPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        thetaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "Theta",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless, 0)
            )
        );
    }
}


template<class CloudType>
void Foam::VoidFractionSC<CloudType>::postEvolve()
{
    volScalarField& theta = thetaPtr_();

    const fvMesh& mesh = this->owner().mesh();

    theta.primitiveFieldRef() /= mesh.time().deltaTValue()*mesh.V();

    CloudFunctionObjectSC<CloudType>::postEvolve();
}


template<class CloudType>
void Foam::VoidFractionSC<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point&,
    bool&
)
{
    volScalarField& theta = thetaPtr_();

    theta[p.cell()] += dt*p.nParticle()*p.volume();
}


// ************************************************************************* //
