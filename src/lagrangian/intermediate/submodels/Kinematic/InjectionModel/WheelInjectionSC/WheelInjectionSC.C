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

#include "WheelInjectionSC.H"
#include "TimeFunction1.H"
#include "Constant.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::WheelInjectionSC<CloudType>::setInjectionMethod()
{
    const word injectionMethod =
        this->coeffDict().template lookupOrDefault<word>
        (
            "injectionMethod",
            word::null
        );

    if (injectionMethod == "tread" || injectionMethod == word::null)
    {
        injectionMethod_ = imTread;
    }
    else if (injectionMethod == "capillary")
    {
        injectionMethod_ = imCapillary;
    }
    else
    {
        FatalErrorInFunction
            << "injectionMethod must be either 'tread' or 'capillary'"
            << exit(FatalError);
    }

    this->coeffDict().lookup("degMin") >> degMin_;
    this->coeffDict().lookup("degMax") >> degMax_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WheelInjectionSC<CloudType>::WheelInjectionSC
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModelSC<CloudType>(dict, owner, modelName, typeName),
    injectionMethod_(imTread),
    position_
    (
        TimeFunction1<vector>
        (
            owner.db().time(),
            "position",
            this->coeffDict()
        )
    ),
    injectorCell_(-1),
    injectorTetFace_(-1),
    injectorTetPt_(-1),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    flowRateProfile_
    (
        TimeFunction1<scalar>
        (
            owner.db().time(),
            "flowRateProfile",
            this->coeffDict()
        )
    ),
    width_(readScalar(this->coeffDict().lookup("width"))),
    radius_(readScalar(this->coeffDict().lookup("radius"))),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"), owner.rndGen()
        )
    ),
    degMin_(0.0),
    degMax_(360.0),
    degInj_(0.0),
    Umag_(owner.db().time(), "Umag")
{
    duration_ = owner.db().time().userTimeToTime(duration_);

    setInjectionMethod();

    Umag_.reset(this->coeffDict());

    // Set total volume to inject
    this->volumeTotal_ = flowRateProfile_.integrate(0, duration_);
}


template<class CloudType>
Foam::WheelInjectionSC<CloudType>::WheelInjectionSC
(
    const WheelInjectionSC<CloudType>& im
)
:
    InjectionModelSC<CloudType>(im),
    injectionMethod_(im.injectionMethod_),
    position_(im.position_),
    injectorCell_(im.injectorCell_),
    injectorTetFace_(im.injectorTetFace_),
    injectorTetPt_(im.injectorTetPt_),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    flowRateProfile_(im.flowRateProfile_),
    width_(im.width_),
    radius_(im.radius_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    degMin_(im.degMin_),
    degMax_(im.degMax_),
    degInj_(im.degInj_),
    Umag_(im.Umag_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WheelInjectionSC<CloudType>::~WheelInjectionSC()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::WheelInjectionSC<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::WheelInjectionSC<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (time0 >= 0 && time0 < duration_)
    {
        //// Standard calculation
        //return floor(parcelsPerSecond_*(time1 - time0));

        // Modified calculation to make numbers exact
        return floor(parcelsPerSecond_*time1 - this->parcelsAddedTotal());
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::WheelInjectionSC<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if (time0 >= 0 && time0 < duration_)
    {
        return flowRateProfile_.integrate(time0, time1);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
void Foam::WheelInjectionSC<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    // Generate random numbers to determine the injection position
    Random& rndGen = this->owner().rndGen();
    scalar fracA = 0;
    scalar fracB = 0;

    if (Pstream::master())
    {
        fracA = rndGen.scalar01();
        fracB = rndGen.scalar01();
    }

    reduce(fracA, maxOp<scalar>());
    reduce(fracB, maxOp<scalar>());

    // Load the ‘wheel’ centre point
    const scalar t = time - this->SOI_;
    position = position_.value(t);

    // Set the injection position in X and Z
    const scalar degDelta = degMax_ - degMin_;
    degInj_ = degToRad(degMin_ + (fracA * degDelta));
    position.component(0) = position.component(0) + (radius_ * cos(degInj_));
    position.component(2) = position.component(2) + (radius_ * sin(degInj_));

    // Set the injection position in Y
    const scalar widthInj = (fracB * width_) - (width_ / 2);
    position.component(1) = position.component(1) + widthInj;

    // Assign the parcel a cell
    this->findCellAtPosition
    (
        cellOwner,
        tetFacei,
        tetPti,
        position,
        false
    );
}


template<class CloudType>
void Foam::WheelInjectionSC<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    const scalar t = time - this->SOI_;

        switch (injectionMethod_)
    {
        case imTread:
        {
            // Calculate the direction of injection
            vector dirVec = vector::max;
            dirVec.component(0) = cos(degInj_ - degToRad(degMin_));
            dirVec.component(1) = 0;
            dirVec.component(2) = sin(degInj_ - degToRad(degMin_));
            dirVec = normalised(dirVec);

            // Set the velocity
            parcel.U() = Umag_.value(t)*dirVec;
            break;
        }
        case imCapillary:
        {
            // Calculate the direction of injection
            vector dirVec = vector::max;
            dirVec.component(0) = -sin(degInj_);
            dirVec.component(1) = 0;
            dirVec.component(2) = cos(degInj_);
            dirVec = normalised(dirVec);

            // Set the velocity
            parcel.U() = Umag_.value(t)*dirVec;
            break;
        }
        default:
        {
            break;
        }
    }

    // Set the particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::WheelInjectionSC<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::WheelInjectionSC<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
