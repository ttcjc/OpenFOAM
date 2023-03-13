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

#include "forceCoeffsExtended.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forceCoeffsExtended, 0);
    addToRunTimeSelectionTable(functionObject, forceCoeffsExtended, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::forceCoeffsExtended::writeFileHeader(const label i)
{
    switch (fileID(i))
    {
        case fileID::mainFile:
        {
            // force coeff data

            writeHeader(file(i), "Force coefficients");
            writeHeaderValue(file(i), "liftDir", liftDir_);
            writeHeaderValue(file(i), "dragDir", dragDir_);
            writeHeaderValue(file(i), "sideDir", sideDir_);
            writeHeaderValue(file(i), "pitchAxis", pitchAxis_);
            writeHeaderValue(file(i), "yawAxis", yawAxis_);
            writeHeaderValue(file(i), "rollAxis", rollAxis_);
            writeHeaderValue(file(i), "magUInf", magUInf_);
            writeHeaderValue(file(i), "lRef", lRef_);
            writeHeaderValue(file(i), "Aref", Aref_);
            writeHeaderValue(file(i), "CofR", coordSys_.origin());
            writeCommented(file(i), "Time");
            writeTabbed(file(i), "Cl");
            writeTabbed(file(i), "Cd");
            writeTabbed(file(i), "Cs");
            writeTabbed(file(i), "Cm(p)");
            writeTabbed(file(i), "Cm(y)");
            writeTabbed(file(i), "Cm(r)");

            break;
        }
        case fileID::binsFile:
        {
            // bin coeff data

            writeHeader(file(i), "Force coefficient bins");
            writeHeaderValue(file(i), "bins", nBin_);
            writeHeaderValue(file(i), "start", binMin_);
            writeHeaderValue(file(i), "delta", binDx_);
            writeHeaderValue(file(i), "direction", binDir_);

            vectorField binPoints(nBin_);
            writeCommented(file(i), "x co-ords  :");
            forAll(binPoints, pointi)
            {
                binPoints[pointi] = (binMin_ + (pointi + 1)*binDx_)*binDir_;
                file(i) << tab << binPoints[pointi].x();
            }
            file(i) << nl;

            writeCommented(file(i), "y co-ords  :");
            forAll(binPoints, pointi)
            {
                file(i) << tab << binPoints[pointi].y();
            }
            file(i) << nl;

            writeCommented(file(i), "z co-ords  :");
            forAll(binPoints, pointi)
            {
                file(i) << tab << binPoints[pointi].z();
            }
            file(i) << nl;

            writeCommented(file(i), "Time");

            for (label j = 0; j < nBin_; j++)
            {
                const word jn('(' + Foam::name(j) + ')');
                writeTabbed(file(i), "Cl" + jn);
                writeTabbed(file(i), "Cd" + jn);
                writeTabbed(file(i), "Cs" + jn);
                writeTabbed(file(i), "Cm(p)" + jn);
                writeTabbed(file(i), "Cm(y)" + jn);
                writeTabbed(file(i), "Cm(r)" + jn);
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }

    file(i)<< endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffsExtended::forceCoeffsExtended
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces(name, runTime, dict),
    liftDir_(0,0,1),
    dragDir_(1,0,0),
    sideDir_(0,1,0),
    pitchAxis_(0,1,0),
    yawAxis_(0,0,1),
    rollAxis_(1,0,0),
    magUInf_(1.0),
    lRef_(1.0),
    Aref_(1.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forceCoeffsExtended::~forceCoeffsExtended()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forceCoeffsExtended::read(const dictionary& dict)
{
    forces::read(dict);

    // Directions for lift, drag and side forces, and pitch, yaw and roll moments
    // Normalize to ensure that the directions are unit vectors
    dict.lookup("liftDir") >> liftDir_;
    liftDir_ /= mag(liftDir_);

    dict.lookup("dragDir") >> dragDir_;
    dragDir_ /= mag(dragDir_);

    dict.lookup("sideDir") >> sideDir_;
    sideDir_ /= mag(sideDir_);

    dict.lookup("pitchAxis") >> pitchAxis_;
    pitchAxis_ /= mag(pitchAxis_);

    dict.lookup("yawAxis") >> yawAxis_;
    yawAxis_ /= mag(yawAxis_);

    dict.lookup("rollAxis") >> rollAxis_;
    rollAxis_ /= mag(rollAxis_);

    // Freestream velocity magnitude
    dict.lookup("magUInf") >> magUInf_;

    // Reference (freestream) density
    dict.lookup("rhoInf") >> rhoRef_;

    // Reference length and area scales
    dict.lookup("lRef") >> lRef_;
    dict.lookup("Aref") >> Aref_;

    return true;
}


bool Foam::functionObjects::forceCoeffsExtended::execute()
{
    return true;
}


bool Foam::functionObjects::forceCoeffsExtended::write()
{
    forces::calcForcesMoment();

    if (Pstream::master())
    {
        logFiles::write();

        scalar pDyn = 0.5*rhoRef_*magUInf_*magUInf_;

        Field<vector> totForce(force_[0] + force_[1] + force_[2]);
        Field<vector> totMoment(moment_[0] + moment_[1] + moment_[2]);

        List<Field<scalar>> coeffs(6);
        coeffs[0].setSize(nBin_);
        coeffs[1].setSize(nBin_);
        coeffs[2].setSize(nBin_);
        coeffs[3].setSize(nBin_);
        coeffs[4].setSize(nBin_);
        coeffs[5].setSize(nBin_);

        // lift, drag and side force and pitch, yaw and roll moments
        coeffs[0] = (totForce & liftDir_)/(Aref_*pDyn);
        coeffs[1] = (totForce & dragDir_)/(Aref_*pDyn);
        coeffs[2] = (totForce & sideDir_)/(Aref_*pDyn);

        coeffs[3] = (totMoment & pitchAxis_)/(Aref_*lRef_*pDyn);
        coeffs[4] = (totMoment & yawAxis_)/(Aref_*lRef_*pDyn);
        coeffs[5] = (totMoment & rollAxis_)/(Aref_*lRef_*pDyn);

        scalar Cl = sum(coeffs[0]);
        scalar Cd = sum(coeffs[1]);
        scalar Cs = sum(coeffs[2]);

        scalar Cmp = sum(coeffs[3]);
        scalar Cmy = sum(coeffs[4]);
        scalar Cmr = sum(coeffs[5]);

        writeTime(file(fileID::mainFile));
        file(fileID::mainFile)
            << tab << Cl << tab  << Cd << tab << Cs
            << tab << Cmp << tab << Cmy << tab << Cmr << endl;

        Log << type() << " " << name() << " write:" << nl
            << "    Cl    = " << Cl << nl
            << "    Cd    = " << Cd << nl
            << "    Cs    = " << Cs << nl
            << "    Cm(p) = " << Cmp << nl
            << "    Cm(y) = " << Cmy << nl
            << "    Cm(r) = " << Cmr << endl;

        if (nBin_ > 1)
        {
            if (binCumulative_)
            {
                for (label i = 1; i < coeffs[0].size(); i++)
                {
                    coeffs[0][i] += coeffs[0][i-1];
                    coeffs[1][i] += coeffs[1][i-1];
                    coeffs[2][i] += coeffs[2][i-1];
                    coeffs[3][i] += coeffs[3][i-1];
                    coeffs[4][i] += coeffs[4][i-1];
                    coeffs[5][i] += coeffs[5][i-1];
                }
            }

            writeTime(file(fileID::binsFile));

            forAll(coeffs[0], i)
            {
                file(fileID::binsFile)
                    << tab << coeffs[5][i]
                    << tab << coeffs[4][i]
                    << tab << coeffs[3][i]
                    << tab << coeffs[2][i]
                    << tab << coeffs[1][i]
                    << tab << coeffs[0][i];
            }

            file(fileID::binsFile) << endl;
        }

        Log << endl;
    }

    return true;
}


// ************************************************************************* //
