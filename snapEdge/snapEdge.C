/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    snapEdge

Description
    snap edges to feature lines

Usage
    Options are:

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dbse.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//  Main program:
int main(int argc, char *argv[])
{
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const word oldInstance = mesh.pointsInstance();
    const bool overwrite = args.optionFound("overwrite");

#   include "readSnapEdgeDict.H"

    dbse model(runTime, mesh, snapEdgeDict);

    model.snapEdges();

    if(!overwrite) 
    {
        runTime++;
    } 
    else 
    {
        mesh.setInstance(oldInstance);
    }

    Info << "Writing points to Time = " << runTime.value() << endl;
    mesh.write();

    return 0;
}


// ************************************************************************* //
