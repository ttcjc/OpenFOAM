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
    Foam::WheelInjectionSC

Description
    A custom injector that introduces particles randomly about a user-defined
    cylinder representing a wheel.

    Two injection methods are available, intended to produce injection patterns
    representative of both tread pickup and capillary action, as defined by
    Weir et al. in "Reduction of Adverse Aerodynamic Effects of Large Trucks"

Usage
    \table
    Property        | Description                    | Required | Default
    injectionMethod | 'Wheel' spray ejection modes   | no       | tread
    position        | 'Wheel' centre point (m m m)   | yes      |
    width_          | 'Wheel' width (m)              | yes      |
    radius          | 'Wheel' radius (m)             | yes      |
    degMin          | Minimum injection angle (deg)  | no       | 0
    degMax          | Maximum injection angle (deg)  | no       | 360
    \endtable

    Example specification:

    \verbatim
    wheelInjector
    {
        type                wheelInjectionSC;

        // Position
        injectionMethod     tread;
        position            (0 0 0);
        width               1;
        radius              1;
        degMin              0;
        degMax              360;

        // Duration
        SOI                 0;
        duration            1;

        // Quantity
        parcelBasisType     mass;
        massTotal           1; // Total mass injected over duration
        parcelsPerSecond    1;
        flowRateProfile     constant 1;

        // Velocity
        Umag                constant 1;

        // Size
        sizeDistribution
        {
            type    fixedValue;

            fixedValueDistribution
            {
                value   1e-6;
            }
        }
    }
    \endverbatim

SourceFiles
    WheelInjectionSC.C

\*---------------------------------------------------------------------------*/

#ifndef WheelInjectionSC_H
#define WheelInjectionSC_H

#include "InjectionModelSC.H"
#include "distributionModel.H"
#include "TimeFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class WheelInjectionSC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WheelInjectionSC
:
    public InjectionModelSC<CloudType>
{
public:

    //- Injection method enumeration
    enum injectionMethod
    {
        imTread,
        imCapillary
    };


private:

    // Private Data

        //- "Tread" or "Capillary" injection method
        injectionMethod injectionMethod_;

        //- Position of the injector
        const TimeFunction1<vector> position_;

        //- Cell label corresponding to the injector position
        label injectorCell_;

        //- Tet-face label corresponding to the injector position
        label injectorTetFace_;

        //- Tet-point label corresponding to the injector position
        label injectorTetPt_;

        //- Injection duration [s]
        scalar duration_;

        //- Number of parcels to introduce per second
        const label parcelsPerSecond_;

        //- Flow rate profile relative to SOI []
        const TimeFunction1<scalar> flowRateProfile_;

        //- Width of injection wheel [m]
        scalar width_;

        //- Radius of injection wheel [m]
        scalar radius_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;

        //- Minimum potential angle of injection [deg]
        scalar degMin_;

        //- Maximum potential angle of injection [deg]
        scalar degMax_;

        //- Actual angle of injection [deg]
        scalar degInj_;

        //- Parcel velocity [m/s]
        TimeFunction1<scalar> Umag_;


    // Private Member Functions

        //- Set the injection type
        void setInjectionMethod();


public:

    //- Runtime type information
    TypeName("wheelInjectionSC");


    // Constructors

        //- Construct from dictionary
        WheelInjectionSC
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        WheelInjectionSC(const WheelInjectionSC<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModelSC<CloudType>> clone() const
        {
            return autoPtr<InjectionModelSC<CloudType>>
            (
                new WheelInjectionSC<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~WheelInjectionSC();


    // Member Functions

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);



        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WheelInjectionSC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
