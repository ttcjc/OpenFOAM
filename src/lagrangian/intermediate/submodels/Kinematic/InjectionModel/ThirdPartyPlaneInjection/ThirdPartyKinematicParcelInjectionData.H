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
    Foam::ThirdPartyKinematicParcelInjectionData

Description
    Container class to provide injection data for kinematic parcels

SourceFiles
    ThirdPartyKinematicParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef ThirdPartyKinematicParcelInjectionData_H
#define ThirdPartyKinematicParcelInjectionData_H

#include "dictionary.H"
#include "vector.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class ThirdPartyKinematicParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const ThirdPartyKinematicParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    ThirdPartyKinematicParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
               Class ThirdPartyKinematicParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class ThirdPartyKinematicParcelInjectionData
{
protected:

    // Parcel properties

        //- Position [m]
        point x_;

        //- Velocity [m/s]
        vector U_;

        //- Diameter [m]
        scalar d_;

        //- Density [kg/m^3]
        scalar rho_;

        //- Mass flow rate [kg/s]
        scalar mDot_;

		// CJC {
		//- Number of particles contained in parcel [dimless]
		scalar particleCount_;
		// } CJC

		// CJC {
		//- Injection time relative to SOI [s]
		scalar injectionTime_;
		// } CJC


public:

    //- Runtime type information
    TypeName("ThirdPartyKinematicParcelInjectionData");

    // Constructors

        //- Null constructor
        ThirdPartyKinematicParcelInjectionData();

        //- Construct from dictionary
        ThirdPartyKinematicParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        ThirdPartyKinematicParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~ThirdPartyKinematicParcelInjectionData();


    // Access

        //- Return const access to the position
        inline const point& x() const;

        //- Return const access to the velocity
        inline const vector& U() const;

        //- Return const access to the diameter
        inline scalar d() const;

        //- Return const access to the density
        inline scalar rho() const;

        //- Return const access to the mass flow rate
        inline scalar mDot() const;

		// CJC {
		//- Return const access to the injector particle count
		inline scalar particleCount() const;
		// } CJC

		// CJC {
		//- Return const access to the injection time
		inline scalar injectionTime() const;
		// } CJC


    // Edit

        //- Return access to the position
        inline point& x();

        //- Return access to the velocity
        inline vector& U();

        //- Return access to the diameter
        inline scalar& d();

        //- Return access to the density
        inline scalar& rho();

        //- Return access to the mass flow rate
        inline scalar& mDot();

		// CJC {
		//- Return access to the injector particle count
		inline scalar& particleCount();
		// } CJC

		// CJC {
		//- Return access to the injection time
		inline scalar& injectionTime();
		// } CJC


    // Friend Operators

        friend bool operator==
        (
            const ThirdPartyKinematicParcelInjectionData& a,
            const ThirdPartyKinematicParcelInjectionData& b
        )
        {
            NotImplemented;

            return false;
        }

        friend bool operator!=
        (
            const ThirdPartyKinematicParcelInjectionData& a,
            const ThirdPartyKinematicParcelInjectionData& b
        )
        {
            NotImplemented;

            return false;
        }


    // I-O

        //- Ostream operator
        friend Ostream& operator<<
        (
            Ostream& os,
            const ThirdPartyKinematicParcelInjectionData& data
        );

        //- Istream operator
        friend Istream& operator>>
        (
            Istream& is,
            ThirdPartyKinematicParcelInjectionData& data
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThirdPartyKinematicParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
