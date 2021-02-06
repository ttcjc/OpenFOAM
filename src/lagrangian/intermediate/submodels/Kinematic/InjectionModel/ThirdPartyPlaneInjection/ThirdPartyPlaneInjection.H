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
    Foam::ThirdPartyPlaneInjection

Description
    Particle injection sources read from look-up table. Each row corresponds to
    an individual parcel. Used alongside particle extraction functionality to
	replicate Third-Party contaminant sources.

    \verbatim
    (
        (x y z) (u v w) d rho mDot particleCount injectionTime // parcel 1
        (x y z) (u v w) d rho mDot particleCount injectionTime // parcel 2
        ...
        (x y z) (u v w) d rho mDot particleCount injectionTime // parcel N
    );
    \endverbatim

    where:
    \plaintable
        x, y, z       | global cartesian co-ordinates [m]
        u, v, w       | global cartesian velocity components [m/s]
        d             | diameter [m]
        rho           | density [kg/m^3]
        mDot          | mass flow rate [kg/s] **unused**
		particleCount | number of particles contained in parcel [dimless]
		injectionTime | injection time relative to SOI [s]
    \endplaintable

SourceFiles
    ThirdPartyPlaneInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ThirdPartyPlaneInjection_H
#define ThirdPartyPlaneInjection_H

#include "InjectionModel.H"
#include "ThirdPartyKinematicParcelInjectionDataIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class ThirdPartyPlaneInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThirdPartyPlaneInjection
:
    public InjectionModel<CloudType>
{
    // Private Data

        //- Name of file containing injector/parcel data
        const word inputFileName_;

        //- Injection duration - common to all injection sources
        scalar duration_;

        //- Number of parcels per injector - common to all injection sources
        const scalar parcelsPerSecond_;

        //- Flag to indicate to randomise injection positions
        bool randomise_;

        //- List of injectors
        ThirdPartyKinematicParcelInjectionDataIOList injectors_;

        //- List of cell labels corresponding to injector positions
        labelList injectorCells_;

        //- List of tetFace labels corresponding to injector positions
        labelList injectorTetFaces_;

        //- List of tetPt labels corresponding to injector positions
        labelList injectorTetPts_;


public:

    //- Runtime type information
    TypeName("ThirdPartyPlaneInjection");


    // Constructors

        //- Construct from dictionary
        ThirdPartyPlaneInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ThirdPartyPlaneInjection
        (
            const ThirdPartyPlaneInjection<CloudType>& im
        );

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new ThirdPartyPlaneInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ThirdPartyPlaneInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

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

			// CJC {
			//- Version of setPositionAndCell for ThirdPartyPlaneInjection
			virtual void setPositionAndCellTP
			(
				const label parcelI,
				const label nParcels,
				const scalar time,
				label& currentInjector,
				vector& position,
				label& cellOwner,
				label& tetFacei,
				label& tetPti
			);
			// } CJC

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

			// CJC {
			//- Version of setProperties for ThirdPartyPlaneInjection
			virtual void setPropertiesTP
			(
				const label currentInjector,
				typename CloudType::parcelType& parcel
			);
			// } CJC

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
    #include "ThirdPartyPlaneInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
