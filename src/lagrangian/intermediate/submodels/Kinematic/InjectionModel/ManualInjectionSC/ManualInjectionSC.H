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
    Foam::ManualInjectionSC

Description
    Manual injection.

    User specifies:
      - Total mass to inject
      - Parcel positions in file \c positionsFile
      - Initial parcel velocity

    Properties:
      - Parcel diameters obtained by distribution model
      - All parcels introduced at SOI

SourceFiles
    ManualInjectionSC.C

\*---------------------------------------------------------------------------*/

#ifndef ManualInjectionSC_H
#define ManualInjectionSC_H

#include "InjectionModelSC.H"
#include "distributionModel.H"
#include "Switch.H"
#include "GlobalIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ManualInjectionSC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ManualInjectionSC
:
    public InjectionModelSC<CloudType>
{
    // Private Data

        //- Name of file containing positions data
        const word positionsFile_;

        //- Field of parcel positions
        GlobalIOField<vector> positions_;

        //- Field of parcel diameters
        scalarList diameters_;

        //- List of cell labels corresponding to injector positions
        labelList injectorCells_;

        //- List of tetFace labels corresponding to injector positions
        labelList injectorTetFaces_;

        //- List of tetPt labels corresponding to injector positions
        labelList injectorTetPts_;

        //- Initial parcel velocity
        const vector U0_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;

        //- Flag to suppress errors if particle injection site is out-of-bounds
        Switch ignoreOutOfBounds_;


public:

    //- Runtime type information
    TypeName("manualInjectionSC");


    // Constructors

        //- Construct from dictionary
        ManualInjectionSC
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ManualInjectionSC(const ManualInjectionSC<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModelSC<CloudType>> clone() const
        {
            return autoPtr<InjectionModelSC<CloudType>>
            (
                new ManualInjectionSC<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ManualInjectionSC();


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
    #include "ManualInjectionSC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
