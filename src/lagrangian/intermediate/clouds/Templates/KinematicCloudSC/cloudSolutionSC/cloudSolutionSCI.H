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

inline const Foam::fvMesh& Foam::cloudSolutionSC::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary& Foam::cloudSolutionSC::dict() const
{
    return dict_;
}


inline const Foam::Switch Foam::cloudSolutionSC::active() const
{
    return active_;
}


inline const Foam::dictionary& Foam::cloudSolutionSC::sourceTermDict() const
{
    return dict_.subDict("sourceTerms");
}


inline const Foam::dictionary& Foam::cloudSolutionSC::interpolationSchemes() const
{
    return dict_.subDict("interpolationSchemes");
}


inline const Foam::dictionary& Foam::cloudSolutionSC::integrationSchemes() const
{
    return dict_.subDict("integrationSchemes");
}


inline const Foam::Switch Foam::cloudSolutionSC::transient() const
{
    return transient_;
}


inline const Foam::Switch Foam::cloudSolutionSC::steadyState() const
{
    return !transient_;
}


// CJC {
    inline Foam::label Foam::cloudSolutionSC::LagrangianWrite() const
    {
        return LagrangianWrite_;
    }


    inline Foam::scalarList Foam::cloudSolutionSC::extractionPlaneList() const
    {
        return extractionPlaneList_;
    }
// } CJC


inline Foam::label Foam::cloudSolutionSC::calcFrequency() const
{
    return calcFrequency_;
}


inline Foam::scalar Foam::cloudSolutionSC::maxCo() const
{
    return maxCo_;
}


inline Foam::label Foam::cloudSolutionSC::iter() const
{
    return iter_;
}


inline Foam::label Foam::cloudSolutionSC::nextIter()
{
    return ++iter_;
}


inline Foam::scalar Foam::cloudSolutionSC::trackTime() const
{
    return trackTime_;
}


inline Foam::Switch& Foam::cloudSolutionSC::coupled()
{
    return coupled_;
}


inline const Foam::Switch Foam::cloudSolutionSC::coupled() const
{
    return coupled_;
}


inline const Foam::Switch Foam::cloudSolutionSC::cellValueSourceCorrection() const
{
    return cellValueSourceCorrection_;
}


inline Foam::scalar Foam::cloudSolutionSC::maxTrackTime() const
{
    return maxTrackTime_;
}


inline const Foam::Switch Foam::cloudSolutionSC::resetSourcesOnStartup() const
{
    return resetSourcesOnStartup_;
}


// ************************************************************************* //
