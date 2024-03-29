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
    Foam::IOPositionSC

Description
    Helper IO class to read and write particle positions

SourceFiles
    IOPositionSC.C

\*---------------------------------------------------------------------------*/

#ifndef IOPositionSC_H
#define IOPositionSC_H

#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IOPositionSC Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IOPositionSC
:
    public regIOobject
{

    // Private Data

        //- Reference to the cloud
        const CloudType& cloud_;


public:

    // Static data

        //- Runtime type name information. Use cloud type.
        virtual const word& type() const
        {
            return CloudSC<typename CloudType::particleType>::typeName;
        }


    // Constructors

        //- Construct from cloud
        IOPositionSC(const CloudType&);


    // Member Functions

        //- Inherit readData from regIOobject
        using regIOobject::readData;

        virtual void readData(Istream&, CloudType&);

        virtual bool write(const bool write = true) const;

        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOPositionSC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
