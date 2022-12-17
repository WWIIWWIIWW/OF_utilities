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

Application
    CalSVR

Description
    Calculate surface to volume ratio and characteristics length (inverse) of structure by kai.zhang.1@city.ac.uk
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    instantList timeDirs = timeSelector::select0(runTime, args);
    
    //OFstream file("time_list.txt");

    //timeDirs not really needed.
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
               
        //Calculate total volume of mesh
        const scalarField cellVol = mesh.V();
        scalar volTot = 0.0;
        
        forAll(cellVol, cellVolI)
        {
            volTot += cellVol[cellVolI];
        }
        
        reduce(volTot, sumOp<scalar>());
        Info<< "Total volume = "
            << volTot
            << endl;

        //Calculate patch total area
        const surfaceScalarField& magSf = mesh.magSf();
        
        const fvPatchList& patches = mesh.boundary();
        
        scalar areaTot = 0.0;
        forAll(patches,patchI) 
        {
            const fvPatch& cPatch = patches[patchI];
            
            scalar patchArea = 0.0;  
            
            forAll(cPatch, faceI)
            {
                patchArea += magSf.boundaryField()[patchI][faceI];
            }
            
            Info << "Patch Area for " << cPatch.name()
                 << "=" << patchArea << endl;
            
            areaTot += patchArea;
        }
        reduce(areaTot, sumOp<scalar>());
        
        Info<< "Total area = "
            << areaTot
            << endl;      
        
        Info<< "SVR = " << areaTot/volTot << endl;
        Info<< "Characteristic Length = " << volTot/areaTot << endl;
           
    }


    Info<< "\nEnd\n" << endl;
    
    return 0;
}


// ************************************************************************* //
