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
    Foam::CloudFunctionObjectSC

Description
    Templated cloud function object base class

SourceFiles
    CloudFunctionObjectSC.C
    CloudFunctionObjectNewSC.C

\*---------------------------------------------------------------------------*/

#ifndef CloudFunctionObjectSC_H
#define CloudFunctionObjectSC_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class tetIndices;

/*---------------------------------------------------------------------------*\
                    Class CloudFunctionObjectSC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CloudFunctionObjectSC
:
    public CloudSubModelBase<CloudType>
{
    // Private Data

        //- Output path
        fileName outputDir_;

        // CJC {
            //- Write interval (match 'LagrangianWrite' in 'kinematicCloudProperties')
            label writeInterval_;
        // } CJC


    // Private Member Functions

        //- Inherite write from CloudSubModelBase
        using CloudSubModelBase<CloudType>::write;

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("CloudFunctionObjectSC");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        CloudFunctionObjectSC,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        ),
        (dict, owner, modelName)
    );


    // Constructors

        //- Construct null from owner
        CloudFunctionObjectSC(CloudType& owner);

        //- Construct from dictionary
        CloudFunctionObjectSC
        (
            const dictionary& dict,
            CloudType& owner,
            const word& objectType,
            const word& modelName
        );

        //- Construct copy
        CloudFunctionObjectSC(const CloudFunctionObjectSC<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObjectSC<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObjectSC<CloudType>>
            (
                new CloudFunctionObjectSC<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CloudFunctionObjectSC();


    //- Selector
    static autoPtr<CloudFunctionObjectSC<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner,
        const word& objectType,
        const word& modelName
    );


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-evolve hook
            virtual void postEvolve();

            //- Post-move hook
            virtual void postMove
            (
                typename CloudType::parcelType& p,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );

            //- Post-patch hook
            virtual void postPatch
            (
                const typename CloudType::parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );

            //- Post-face hook
            virtual void postFace
            (
                const typename CloudType::parcelType& p,
                bool& keepParticle
            );


        // Input/output

            //- Return the output path
            const fileName& outputDir() const;

            //- Return the output time path
            fileName writeTimeDir() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCloudFunctionObjectSC(CloudType)                                   \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::CloudFunctionObjectSC<kinematicCloudType>,                       \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            CloudFunctionObjectSC<kinematicCloudType>,                         \
            dictionary                                                         \
        );                                                                     \
    }


#define makeCloudFunctionObjectSCType(SS, CloudType)                           \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<kinematicCloudType>, 0);      \
                                                                               \
    Foam::CloudFunctionObjectSC<kinematicCloudType>::                          \
        adddictionaryConstructorToTable<Foam::SS<kinematicCloudType>>          \
            add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CloudFunctionObjectSC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
