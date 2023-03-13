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

#include "CloudFunctionObjectSC.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudFunctionObjectSC<CloudType>::write()
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudFunctionObjectSC<CloudType>::CloudFunctionObjectSC(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    outputDir_()
{}


template<class CloudType>
Foam::CloudFunctionObjectSC<CloudType>::CloudFunctionObjectSC
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName,
    const word& objectType
)
:
    CloudSubModelBase<CloudType>(modelName, owner, dict, typeName, objectType),
    outputDir_(owner.mesh().time().path())
{
    const fileName relPath =
        "postProcessing"/cloud::prefix/owner.name()/this->modelName();


    if (Pstream::parRun())
    {
        // Put in undecomposed case (Note: gives problems for
        // distributed data running)
        outputDir_ = outputDir_/".."/relPath;
    }
    else
    {
        outputDir_ = outputDir_/relPath;
    }
    outputDir_.clean();

    // CJC {
        // Store Cloud Function Object write interval
        dict.readIfPresent("writeInterval", writeInterval_);
    // } CJC

}


template<class CloudType>
Foam::CloudFunctionObjectSC<CloudType>::CloudFunctionObjectSC
(
    const CloudFunctionObjectSC<CloudType>& ppm
)
:
    CloudSubModelBase<CloudType>(ppm),
    outputDir_(ppm.outputDir_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CloudFunctionObjectSC<CloudType>::~CloudFunctionObjectSC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CloudFunctionObjectSC<CloudType>::preEvolve()
{}


template<class CloudType>
void Foam::CloudFunctionObjectSC<CloudType>::postEvolve()
{
    // CJC {
        // Write Cloud Function Object
        if ((this->owner().time().writeTime()) || (this->owner().time().timeIndex() % writeInterval_ == 0))
        {
            this->write();
        }
    // } CJC
}


template<class CloudType>
void Foam::CloudFunctionObjectSC<CloudType>::postMove
(
    typename CloudType::parcelType&,
    const scalar,
    const point&,
    bool&
)
{}


template<class CloudType>
void Foam::CloudFunctionObjectSC<CloudType>::postPatch
(
    const typename CloudType::parcelType&,
    const polyPatch&,
    bool&
)
{}


template<class CloudType>
void Foam::CloudFunctionObjectSC<CloudType>::postFace
(
    const typename CloudType::parcelType&,
    bool&
)
{}


template<class CloudType>
const Foam::fileName& Foam::CloudFunctionObjectSC<CloudType>::outputDir() const
{
    return outputDir_;
}


template<class CloudType>
Foam::fileName Foam::CloudFunctionObjectSC<CloudType>::writeTimeDir() const
{
    return outputDir_/this->owner().time().timeName();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CloudFunctionObjectNewSC.C"

// ************************************************************************* //
