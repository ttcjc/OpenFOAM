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
    Foam::ParticleTracksSC

Description
    Records particle state (all variables) on each call to postFace

SourceFiles
    ParticleTracksSC.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleTracksSC_H
#define ParticleTracksSC_H

#include "CloudFunctionObjectSC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleTracksSC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleTracksSC
:
    public CloudFunctionObjectSC<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;

            //- Convenience typedef for hash hit-table
            typedef HashTable<label, labelPair, typename labelPair::Hash<>>
                hitTableType;

        //- Number of face-hit intervals between storing parcel data
        label trackInterval_;

        //- Maximum number of particles to store per track
        label maxSamples_;

        //- Flag to indicate whether data should be reset/cleared on writing
        Switch resetOnWrite_;

        //- Table of number of times a particle has hit a face
        hitTableType faceHitCounter_;

        //- Pointer to the cloud storage
        autoPtr<CloudSC<parcelType>> cloudPtr_;


protected:

    // Protected member functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("particleTracksSC");


    // Constructors

        //- Construct from dictionary
        ParticleTracksSC
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleTracksSC(const ParticleTracksSC<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObjectSC<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObjectSC<CloudType>>
            (
                new ParticleTracksSC<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleTracksSC();


    // Member Functions

        // Access

            //- Return const access to the track interval
            inline label trackInterval() const;

            //- Return const access to the max samples
            inline label maxSamples() const;

            //- Return const access to the reset on write flag
            inline const Switch& resetOnWrite() const;

            //- Return the table of number of times a particle has hit a face
            inline const hitTableType& faceHitCounter() const;

            //- Return const access to the cloud
            inline const CloudSC<parcelType>& cloud() const;


        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-face hook
            virtual void postFace(const parcelType& p, bool& keepParticle);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleTracksSCI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleTracksSC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
