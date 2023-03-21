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

#ifndef makeParcelInjectionModelsSC_H
#define makeParcelInjectionModelsSC_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CellZoneInjectionSC.H"
#include "ConeInjectionSC.H"
#include "FieldActivatedInjectionSC.H"
#include "InflationInjectionSC.H"
#include "KinematicLookupTableInjectionSC.H"
#include "ManualInjectionSC.H"
#include "NoInjectionSC.H"
#include "PatchInjectionSC.H"
#include "PatchFlowRateInjectionSC.H"
#include "WheelInjectionSC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelInjectionModelsSC(CloudType)                                 \
                                                                               \
    makeInjectionModelSC(CloudType);                                           \
                                                                               \
    makeInjectionModelSCType(CellZoneInjectionSC, CloudType);                  \
    makeInjectionModelSCType(ConeInjectionSC, CloudType);                      \
    makeInjectionModelSCType(FieldActivatedInjectionSC, CloudType);            \
    makeInjectionModelSCType(InflationInjectionSC, CloudType);                 \
    makeInjectionModelSCType(KinematicLookupTableInjectionSC, CloudType);      \
    makeInjectionModelSCType(ManualInjectionSC, CloudType);                    \
    makeInjectionModelSCType(NoInjectionSC, CloudType);                        \
    makeInjectionModelSCType(PatchInjectionSC, CloudType);                     \
    makeInjectionModelSCType(PatchFlowRateInjectionSC, CloudType);             \
    makeInjectionModelSCType(WheelInjectionSC, CloudType);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
