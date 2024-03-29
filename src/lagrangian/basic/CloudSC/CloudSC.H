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
    Foam::CloudSC

Description
    Base cloud calls templated on particle type

SourceFiles
    CloudSC.C
    CloudSCIO.C

\*---------------------------------------------------------------------------*/

#ifndef CloudSC_H
#define CloudSC_H

#include "cloud.H"
#include "IDLList.H"
#include "IOField.H"
#include "CompactIOField.H"
#include "polyMesh.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of functions
template<class ParticleType>
class CloudSC;

template<class ParticleType>
class IOPositionSC;

template<class ParticleType>
Ostream& operator<<
(
    Ostream&,
    const CloudSC<ParticleType>&
);


/*---------------------------------------------------------------------------*\
                           Class CloudSC Declaration
\*---------------------------------------------------------------------------*/

template<class ParticleType>
class CloudSC
:
    public cloud,
    public IDLList<ParticleType>
{
    // Private Data

        //- Reference to the mesh
        const polyMesh& polyMesh_;

        //- Temporary storage for the global particle positions
        mutable autoPtr<vectorField> globalPositionsPtr_;


    // Private Member Functions

        //- Check patches
        void checkPatches() const;

        //- Initialise cloud on IO constructor
        void initCloud(const bool checkClass);

        //- Read cloud properties dictionary
        void readCloudUniformProperties();

        //- Write cloud properties dictionary
        void writeCloudUniformProperties() const;


public:

    friend class particle;
    template<class ParticleT>
    friend class IOPositionSC;

    typedef ParticleType particleType;

    typedef typename IDLList<ParticleType>::iterator iterator;
    typedef typename IDLList<ParticleType>::const_iterator const_iterator;

    //-Runtime type information
    TypeName("CloudSC");


    // Static data

        //- Name of cloud properties dictionary
        static word cloudPropertiesName;


    // Constructors

        //- Construct from mesh and a list of particles
        CloudSC
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<ParticleType>& particles
        );

        //- Construct from mesh by reading from file with given cloud instance
        //  Optionally disable checking of class name for post-processing
        CloudSC
        (
            const polyMesh& pMesh,
            const word& cloudName,
            const bool checkClass = true
        );


    // Member Functions

        // Access

            //- Return the polyMesh reference
            const polyMesh& pMesh() const
            {
                return polyMesh_;
            }

            //- Return the number of particles in the cloud
            label size() const
            {
                return IDLList<ParticleType>::size();
            };


            // Iterators

                const const_iterator begin() const
                {
                    return IDLList<ParticleType>::begin();
                };

                const const_iterator cbegin() const
                {
                    return IDLList<ParticleType>::cbegin();
                };

                const const_iterator end() const
                {
                    return IDLList<ParticleType>::end();
                };

                const const_iterator cend() const
                {
                    return IDLList<ParticleType>::cend();
                };

                iterator begin()
                {
                    return IDLList<ParticleType>::begin();
                };

                iterator end()
                {
                    return IDLList<ParticleType>::end();
                };


        // Edit

            void clear()
            {
                IDLList<ParticleType>::clear();
            };

            //- Transfer particle to cloud
            void addParticle(ParticleType* pPtr);

            //- Remove particle from cloud and delete
            void deleteParticle(ParticleType&);

            //- Remove lost particles from cloud and delete
            void deleteLostParticles();

            //- Reset the particles
            void cloudReset(const CloudSC<ParticleType>& c);

            //- Move the particles
            template<class TrackCloudType>
            void move
            (
                TrackCloudType& cloud,
                typename ParticleType::trackingData& td,
                const scalar trackTime,
                const scalarList extractionPlanePosition,
                const scalarList geometryBoundingBox
            );

            //- Remap the cells of particles corresponding to the
            //  mesh topology change
            void autoMap(const mapPolyMesh&);


        // Read

            //- Helper to construct IOobject for field and current time.
            IOobject fieldIOobject
            (
                const word& fieldName,
                const IOobject::readOption r
            ) const;

            //- Check lagrangian data field
            template<class DataType>
            void checkFieldIOobject
            (
                const CloudSC<ParticleType>& c,
                const IOField<DataType>& data
            ) const;

            //- Check lagrangian data fieldfield
            template<class DataType>
            void checkFieldFieldIOobject
            (
                const CloudSC<ParticleType>& c,
                const CompactIOField<Field<DataType>, DataType>& data
            ) const;


        // Write

            //- Write the field data for the cloud of particles Dummy at
            //  this level.
            virtual void writeFields() const;

            //- Write using given format, version and compression.
            //  Only writes the cloud file if the Cloud isn't empty
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write = true
            ) const;

            //- Write positions to \<cloudName\>_positions.obj file
            void writePositions() const;

            //- Call this before a topology change. Stores the particles global
            //  positions in the database for use during mapping.
            void storeGlobalPositions() const;

            // CJC {
                // Particle Data Extraction

                    //- Write extraction plane data
                    void extractPlaneData(ParticleType& p, scalar planePos);

                    //- Write surface contamination data
                    void extractSoilingData(ParticleType& p);
            // } CJC


    // Ostream Operator

        friend Ostream& operator<< <ParticleType>
        (
            Ostream&,
            const CloudSC<ParticleType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudSC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
