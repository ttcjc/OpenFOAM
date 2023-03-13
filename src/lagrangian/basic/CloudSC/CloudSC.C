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

\*---------------------------------------------------------------------------*/

#include "CloudSC.H"
#include "processorPolyPatch.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "mapPolyMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "wallPolyPatch.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::CloudSC<ParticleType>::checkPatches() const
{
    const polyBoundaryMesh& pbm = polyMesh_.boundaryMesh();
    bool ok = true;
    forAll(pbm, patchi)
    {
        if (isA<cyclicAMIPolyPatch>(pbm[patchi]))
        {
            const cyclicAMIPolyPatch& cami =
                refCast<const cyclicAMIPolyPatch>(pbm[patchi]);

            ok = ok && cami.singlePatchProc() != -1;
        }
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "Particle tracking across AMI patches is only currently "
            << "supported for cases where the AMI patches reside on a "
            << "single processor" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::CloudSC<ParticleType>::CloudSC
(
    const polyMesh& pMesh,
    const word& cloudName,
    const IDLList<ParticleType>& particles
)
:
    cloud(pMesh, cloudName),
    IDLList<ParticleType>(),
    polyMesh_(pMesh),
    globalPositionsPtr_()
{
    checkPatches();

    // Ask for the tetBasePtIs and oldCellCentres to trigger all processors to
    // build them, otherwise, if some processors have no particles then there
    // is a comms mismatch.
    polyMesh_.tetBasePtIs();
    polyMesh_.oldCellCentres();

    if (particles.size())
    {
        IDLList<ParticleType>::operator=(particles);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
void Foam::CloudSC<ParticleType>::addParticle(ParticleType* pPtr)
{
    this->append(pPtr);
}


template<class ParticleType>
void Foam::CloudSC<ParticleType>::deleteParticle(ParticleType& p)
{
    delete(this->remove(&p));
}


template<class ParticleType>
void Foam::CloudSC<ParticleType>::deleteLostParticles()
{
    forAllIter(typename CloudSC<ParticleType>, *this, pIter)
    {
        ParticleType& p = pIter();

        if (p.cell() == -1)
        {
            WarningInFunction
                << "deleting lost particle at position " << p.position()
                << endl;

            deleteParticle(p);
        }
    }
}


template<class ParticleType>
void Foam::CloudSC<ParticleType>::cloudReset(const CloudSC<ParticleType>& c)
{
    // Reset particle count and particles only
    // - not changing the cloud object registry or reference to the polyMesh
    ParticleType::particleCount_ = 0;
    IDLList<ParticleType>::operator=(c);
}


template<class ParticleType>
template<class TrackCloudType>
void Foam::CloudSC<ParticleType>::move
(
    TrackCloudType& cloud,
    typename ParticleType::trackingData& td,
    const scalar trackTime,
    const scalarList extractionPlanePosition,
    const scalarList geometryBoundingBox
)
{
    const polyBoundaryMesh& pbm = pMesh().boundaryMesh();
    const globalMeshData& pData = polyMesh_.globalData();

    // Which patches are processor patches
    const labelList& procPatches = pData.processorPatches();

    // Indexing of equivalent patch on neighbour processor into the
    // procPatches list on the neighbour
    const labelList& procPatchNeighbours = pData.processorPatchNeighbours();

    // Which processors this processor is connected to
    const labelList& neighbourProcs = pData[Pstream::myProcNo()];

    // Indexing from the processor number into the neighbourProcs list
    labelList neighbourProcIndices(Pstream::nProcs(), -1);

    forAll(neighbourProcs, i)
    {
        neighbourProcIndices[neighbourProcs[i]] = i;
    }

    // Initialise the stepFraction moved for the particles
    forAllIter(typename CloudSC<ParticleType>, *this, pIter)
    {
        pIter().reset();
    }

    // List of lists of particles to be transferred for all of the
    // neighbour processors
    List<IDLList<ParticleType>> particleTransferLists
    (
        neighbourProcs.size()
    );

    // List of destination processorPatches indices for all of the
    // neighbour processors
    List<DynamicList<label>> patchIndexTransferLists
    (
        neighbourProcs.size()
    );

    // Allocate transfer buffers
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Clear the global positions as there are about to change
    globalPositionsPtr_.clear();

    // While there are particles to transfer
    while (true)
    {
        particleTransferLists = IDLList<ParticleType>();
        forAll(patchIndexTransferLists, i)
        {
            patchIndexTransferLists[i].clear();
        }

        // Loop over all particles
        forAllIter(typename CloudSC<ParticleType>, *this, pIter)
        {
            ParticleType& p = pIter();

            // CJC {
                // Store particle position prior to evolution (used for planar extraction)
                vector initialPosition = p.position();
            // } CJC

            // Move the particle
            bool keepParticle = p.move(cloud, td, trackTime);

            // CJC {
                // Check for presence of extraction planes
                if (extractionPlanePosition.size() != 0)
                {
                    // Loop over extraction planes
                    forAll(extractionPlanePosition, i)
                    {
                        // Check if particle has crossed current extraction plane
                        if
                        (
                            (p.active() && (initialPosition.x() < extractionPlanePosition[i] && p.position().x() > extractionPlanePosition[i])) ||
                            (p.active() && (initialPosition.x() > extractionPlanePosition[i] && p.position().x() < extractionPlanePosition[i]))
                        )
                        {
                            // Extract particle data from extraction plane
                            extractPlaneData(p, extractionPlanePosition[i]);

                            // Check if particle has crossed final extraction plane
                            if (i == extractionPlanePosition.size() - 1)
                            {
                                // Flag particle for deletion
                                keepParticle = false;
                            }
                        }
                    }
                }

                // Check for presence of valid geometry bounding box
                if (geometryBoundingBox.size() == 6)
                {
                    // Check if particle is inactive (i.e. it has struck the test geometry)
                    if
                    (
                        !p.active() &&
                        (p.position().x() > (geometryBoundingBox[0] - 1e-7) && p.position().x() < (geometryBoundingBox[1] + 1e-7)) &&
                        (p.position().y() > (geometryBoundingBox[2] - 1e-7) && p.position().y() < (geometryBoundingBox[3] + 1e-7)) &&
                        (p.position().z() > (geometryBoundingBox[4] - 1e-7) && p.position().z() < (geometryBoundingBox[5] + 1e-7))
                    )
                    {
                        // Extract particle data from test geometry
                        extractSoilingData(p);

                        // Flag particle for deletion
                        keepParticle = false;
                    }
                }
            // } CJC

            // If the particle is to be kept
            // (i.e. it hasn't passed through an inlet or outlet)
            if (keepParticle)
            {
                if (td.switchProcessor)
                {
                    #ifdef FULLDEBUG
                    if
                    (
                        !Pstream::parRun()
                     || !p.onBoundaryFace()
                     || procPatchNeighbours[p.patch()] < 0
                    )
                    {
                        FatalErrorInFunction
                            << "Switch processor flag is true when no parallel "
                            << "transfer is possible. This is a bug."
                            << exit(FatalError);
                    }
                    #endif

                    const label patchi = p.patch();

                    const label n = neighbourProcIndices
                    [
                        refCast<const processorPolyPatch>
                        (
                            pbm[patchi]
                        ).neighbProcNo()
                    ];

                    p.prepareForParallelTransfer();

                    particleTransferLists[n].append(this->remove(&p));

                    patchIndexTransferLists[n].append
                    (
                        procPatchNeighbours[patchi]
                    );
                }
            }
            else
            {
                deleteParticle(p);
            }
        }

        if (!Pstream::parRun())
        {
            break;
        }


        // Clear transfer buffers
        pBufs.clear();

        // Stream into send buffers
        forAll(particleTransferLists, i)
        {
            if (particleTransferLists[i].size())
            {
                UOPstream particleStream
                (
                    neighbourProcs[i],
                    pBufs
                );

                particleStream
                    << patchIndexTransferLists[i]
                    << particleTransferLists[i];
            }
        }


        // Start sending. Sets number of bytes transferred
        labelList allNTrans(Pstream::nProcs());
        pBufs.finishedSends(allNTrans);


        bool transferred = false;

        forAll(allNTrans, i)
        {
            if (allNTrans[i])
            {
                transferred = true;
                break;
            }
        }
        reduce(transferred, orOp<bool>());

        if (!transferred)
        {
            break;
        }

        // Retrieve from receive buffers
        forAll(neighbourProcs, i)
        {
            label neighbProci = neighbourProcs[i];

            label nRec = allNTrans[neighbProci];

            if (nRec)
            {
                UIPstream particleStream(neighbProci, pBufs);

                labelList receivePatchIndex(particleStream);

                IDLList<ParticleType> newParticles
                (
                    particleStream,
                    typename ParticleType::iNew(polyMesh_)
                );

                label pI = 0;

                forAllIter(typename CloudSC<ParticleType>, newParticles, newpIter)
                {
                    ParticleType& newp = newpIter();

                    label patchi = procPatches[receivePatchIndex[pI++]];

                    newp.correctAfterParallelTransfer(patchi, td);

                    addParticle(newParticles.remove(&newp));
                }
            }
        }
    }
}


template<class ParticleType>
void Foam::CloudSC<ParticleType>::autoMap(const mapPolyMesh& mapper)
{
    if (!globalPositionsPtr_.valid())
    {
        FatalErrorInFunction
            << "Global positions are not available. "
            << "CloudSC::storeGlobalPositions has not been called."
            << exit(FatalError);
    }

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();
    polyMesh_.oldCellCentres();

    const vectorField& positions = globalPositionsPtr_();

    label i = 0;
    forAllIter(typename CloudSC<ParticleType>, *this, iter)
    {
        iter().autoMap(positions[i], mapper);
        ++ i;
    }
}


template<class ParticleType>
void Foam::CloudSC<ParticleType>::writePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/this->name() + "_positions.obj"
    );

    forAllConstIter(typename CloudSC<ParticleType>, *this, pIter)
    {
        const ParticleType& p = pIter();
        pObj<< "v " << p.position().x() << " " << p.position().y() << " "
            << p.position().z() << nl;
    }

    pObj.flush();
}


template<class ParticleType>
void Foam::CloudSC<ParticleType>::storeGlobalPositions() const
{
    // Store the global positions for later use by autoMap. It would be
    // preferable not to need this. If the mapPolyMesh object passed to autoMap
    // had a copy of the old mesh then the global positions could be recovered
    // within autoMap, and this pre-processing would not be necessary.

    globalPositionsPtr_.reset(new vectorField(this->size()));

    vectorField& positions = globalPositionsPtr_();

    label i = 0;
    forAllConstIter(typename CloudSC<ParticleType>, *this, iter)
    {
        positions[i] = iter().position();
        ++ i;
    }
}


// CJC {
    template<class ParticleType>
    void Foam::CloudSC<ParticleType>::extractPlaneData
    (
        ParticleType& p,
        scalar planePos
    )
    {
        word fileName = Foam::name(planePos);

        OFstream planeData
        (
            "LagrangianExtractionPlane/LagrangianExtractionPlaneData_" + fileName,
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true
        );

        planeData << time().value() << " "
                  << p.origId() << " "
                  << p.origProc() << " "
                  << p.d() << " "
                  << p.nParticle() << " "
                  << p.positionCartesian().x() << " "
                  << p.positionCartesian().y() << " "
                  << p.positionCartesian().z() << " "
                  << p.U().x() << " "
                  << p.U().y() << " "
                  << p.U().z() << " "
                  << p.Uslip().x() << " "
                  << p.Uslip().y() << " "
                  << p.Uslip().z() << nl;
    }


    template<class ParticleType>
    void Foam::CloudSC<ParticleType>::extractSoilingData
    (
        ParticleType& p
    )
    {
        OFstream surfaceData
        (
            "LagrangianSurfaceContamination/LagrangianSurfaceContaminationData",
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true
        );

        surfaceData << time().value() << " "
                    << p.origId() << " "
                    << p.origProc() << " "
                    << p.d() << " "
                    << p.nParticle() << " "
                    << p.positionCartesian().x() << " "
                    << p.positionCartesian().y() << " "
                    << p.positionCartesian().z() << " "
                    << p.U().x() << " "
                    << p.U().y() << " "
                    << p.U().z() << " "
                    << p.Uslip().x() << " "
                    << p.Uslip().y() << " "
                    << p.Uslip().z() << nl;
    }
// } CJC


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CloudSCIO.C"

// ************************************************************************* //
