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

#include "kinematicParcelInjectionSCData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::point& Foam::kinematicParcelInjectionSCData::x() const
{
    return x_;
}


inline const Foam::vector& Foam::kinematicParcelInjectionSCData::U() const
{
    return U_;
}


inline Foam::scalar Foam::kinematicParcelInjectionSCData::d() const
{
    return d_;
}


inline Foam::scalar Foam::kinematicParcelInjectionSCData::rho() const
{
    return rho_;
}


inline Foam::scalar Foam::kinematicParcelInjectionSCData::mDot() const
{
    return mDot_;
}


inline Foam::point& Foam::kinematicParcelInjectionSCData::x()
{
    return x_;
}


inline Foam::vector& Foam::kinematicParcelInjectionSCData::U()
{
    return U_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionSCData::d()
{
    return d_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionSCData::rho()
{
    return rho_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionSCData::mDot()
{
    return mDot_;
}


// ************************************************************************* //
