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
    fvc

Description
    $ wmake
    $ fvc > log
    Finite volume method test code with explanation/comments by Kai Zhang

\*---------------------------------------------------------------------------*/

#include "fvCFD.H" //#include several import headers and --> using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //Let's list 
    //or we can use setRootCaseLists.H which is base class declaration
    #include "setRootCase.H" 
    ////Foam::Info<< "Create time\n" << Foam::endl; 
    //Foam::Time runTime(Foam::Time::controlDictName, args);
    #include "createTime.H" 
    // Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr); #create a ptr (Must read)
    // Foam::Info << " for time = " << runTime.timeName() << Foam::nl;
    // Foam::fvMesh& mesh = meshPtr(); #Now we can access mesh which take address of pointer
    #include "createMesh.H"

    //definition volFieldsFwd.H
    
    //typedef GeometricField<scalar, fvPatchField, volMesh> volScalarField;
    //(1)fvPatchField, can be found at fvPatchField.H, templated class using <Type>
    //<Type> can be found in fieldTypes.H (e.g., labelField, scalarField, vectorField, etc.)
    
    //(2) volMesh is constructed from GeoMesh and from meshPtr
    volScalarField P(pow(mesh.C().component(vector::X), 1)); //get mesh centroid vector X, assign value 1
    P.write();
    //Using "fvc::" will calculate the derivative using the current values. 
    //While "fvm::" will discretize the term into the matrix equation.
    //So "fvc::" will generate a field, while "fvm::" will return matrix coefficients. 
    //So, if you're solving [M] [Q] = [B]
    //fvm:: generates coefficients for [M]
    //fvc:: goes in [B]

    volScalarField gradP(fvc::grad(P)().component(vector::X)); //get the X component of gradP
    gradP.write();

    volVectorField curlC(fvc::curl(1.0*mesh.C()));
    curlC.write();

    surfaceScalarField xf(mesh.Cf().component(vector::X)); //get mesh face centroid vector X
    surfaceScalarField xf4(pow(xf, 4));
    
    Info<< "Here is the list containing vortices of face (along X only):" << endl;
    Info<< xf << endl;
    
    Info<< "Give internal and boundary coordinate:" << endl;
    Info<< mesh.C() << endl;
    Info<< "end" << endl;
}


// ************************************************************************* //
