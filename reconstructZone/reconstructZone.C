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

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    instantList timeDirs = timeSelector::select0(runTime, args);
    
    word cellSetName = "c0";

    label zoneID = mesh.cellZones().findZoneID(cellSetName);

    if (zoneID == -1)
    {
        FatalErrorIn("yourFunctionName")
            << "Cannot find cellZone " << cellSetName << endl
            << "Valid cellZones are " << mesh.cellZones().names()
            << exit(FatalError);
    }

    const labelList& cells = mesh.cellZones()[zoneID];

    Info << "Cells in cellzone " << cellSetName << ":" << endl;
    forAll(cells, i)
    {
        const label cell = cells[i];
        Info << cell << endl;
    }



    Info<< "\nEnd\n" << endl;
    
    return 0;
}



