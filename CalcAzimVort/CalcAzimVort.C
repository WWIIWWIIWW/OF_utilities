/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    solutionGradient

Description
    Calculates the gradient for field T in each timestep

\*---------------------------------------------------------------------------*/

#   include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
	#   include "setRootCase.H"
	#   include "createTime.H"
	#   include "createMesh.H"
	#   include "CalcAzimVort.H"


    Info<< "\nCalculating solution gradient\n" << endl;

	instantList timeDirs = timeSelector::select0(runTime, args);	

	forAll(timeDirs, timeI)
    {
		runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

		Info<< "Reading field U\n" << endl;
		

    	
        forAll (mesh.C(), celli)
        {
            radialUnitVector[celli] = mesh.C()[celli];
            radialUnitVector[celli].component(vector::X) = 0; //set x = 0
    	    radialUnitVector[celli] /= mag(radialUnitVector[celli]); //calculated radialUnitVector
        }    	
        Info <<"1"<< endl;
    	volScalarField Ux = U.component(vector::X);
    	volScalarField Uy = U.component(vector::Y);
    	volScalarField Uz = U.component(vector::Z);
    	Info <<"2"<< endl;
    	tmp<volVectorField> tgradUx(fvc::grad(Ux));
    	Info <<"3"<< endl;
        const volVectorField& gradUx = tgradUx();
        Info <<"3.1"<< endl;
        // dux/dr, duy/dr, duz/dr
        tmp<volScalarField> tdUxdr(radialUnitVector & gradUx);   //map gradU along radial direction
        const volScalarField& dUxdr =tdUxdr();
        Info <<"3.2"<< endl;
    	tmp<volScalarField> tUr(Foam::sqrt(magSqr(Uy) + magSqr(Uz))); //sqrt(Uy^2 + Uz^2)
        const volScalarField& Ur = tUr();
        Info <<"4"<< endl;
    	tmp<volVectorField> tgradUr(fvc::grad(Ur)); //dUr/dx, dUr/dy, dUr/dz
        const volVectorField& gradUr = tgradUr();
        Info <<"5"<< endl;
    	tmp<volScalarField> tdUrdx(gradUr.component(vector::X)); //dUr/dx
        const volScalarField& dUrdx = tdUrdx();
        Info <<"6"<< endl;
        AzimVort = dUrdx -  dUxdr;
		AzimVort.write();
	}

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
