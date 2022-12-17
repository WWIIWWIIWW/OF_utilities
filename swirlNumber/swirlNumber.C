/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "swirlNumber.H"
#include "addToRunTimeSelectionTable.H" 
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(swirlNumber, 0);
    addToRunTimeSelectionTable(functionObject, swirlNumber, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
// NOTE: this returns a list of file names which match indices of the enum
// defined in the header of this class. These names are used to create matching
// files by the logFiles object.
Foam::wordList Foam::functionObjects::swirlNumber::createFileNames
(
    const dictionary& dict
) const
{
    DynamicList<word> names(1);

    // use type of the utility as specified in the dict as the top-level dir name
    const word objectType(dict.lookup("type"));

    // Name for file(MAIN_FILE=0)
    names.append(objectType);

    return names;
}

// NOTE: this method first gets declared in logFiles.H, from which this
// class is derived. This method gets called automatically when the base object's
// write() function gets called too.
// The purpose of the function is to add the header to the output data file.
void Foam::functionObjects::swirlNumber::writeFileHeader(const label i)
{
    // Find the correct file to write to from the enum defined in the header.
    switch (fileID(i))
    {
        case MAIN_FILE:
        {
            writeHeader(file(i), "Swirl number across face zone");
            writeHeaderValue(file(i), "Face zone name", faceZoneName_);
            file() << endl;
            break; // exit the case structure
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::swirlNumber::swirlNumber
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    // NOTE: call the base class constructor
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),

    name_(name),
    active_(true),
    velocityFieldName_("U"),
    // NOTE: Read the face zone to integrate over. Get its name from the dict, find
    // it in the mesh, and get a reference to the list of its faces.
    faceZoneName_(dict.lookup("faceZoneName")), 
    // pointZoneName_(dict.lookup("pointZoneName")),
    
    faceZoneLabel_( mesh_.faceZones().findZoneID(faceZoneName_) ),
    // pointZoneLabel_( mesh_.pointZones().findZoneID(pointZoneName_) ),
    
    patchRadius_(0.),
    p_inf(0.), // Rahul:8/2/22: added p_inf    
    patchNormalVector_(vector::zero),
    originVector_(vector::zero), 
   
    faces_(mesh_.faceZones()[faceZoneLabel_])

    // points_(mesh_.pointZones()[pointZoneLabel_])
{
    // NOTE: calls the separate .read() method to import the controls from the dict.
    // dict reference is passed automatically by the OpenFOAM runtime object manager.
    read(dict);

    // built-in logFiles method for creating file streams.
    resetNames(createFileNames(dict));

    if (active_)
    {
        // Any extra initialisation goes here, if necessary

        // NOTE: type() returns the typeName, as defined in the header file. Name
        // is the individual identifier of this instance, as specified in the dict
        Info << "Finished initialising " << type() << ": " << name_ << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::swirlNumber::~swirlNumber()
{}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::functionObjects::swirlNumber::read(const dictionary& dict)
{
    if (active_)
    {
        velocityFieldName_ = dict.lookupOrDefault<word>("velocityFieldName", "U");
	    patchRadius_= readScalar(dict.lookup("radius"));
        p_inf= readScalar(dict.lookup("p_inf")); // Rahul: 8/2/22: Adding p_inf	
        originVector_ = dict.lookup("origin");
        patchNormalVector_ = dict.lookup("normal");
    }
    return true;
}

bool Foam::functionObjects::swirlNumber::execute()
{
    return true;
}

bool Foam::functionObjects::swirlNumber::end()
{
    if (active_)
    {
        execute();
    }
    return true;
}

void Foam::functionObjects::swirlNumber::timeSet()
{}

Foam::vectorField Foam::functionObjects::swirlNumber::faceCentres
(
     const labelList& zfs,
     vectorField& p
)
{
    const label nFaces = zfs.size();
    vectorField fCtrs(nFaces, Zero);
    //surfaceVectorField fCtrs(nFaces, Zero);

    // for(label faceI = 0; faceI != zfs.size(); faceI++)
    forAll(fCtrs, faceI)
    {
	    fCtrs[faceI] = p[zfs[faceI]]; // Rahul: 25/10/21: to calculate the faceCentre.  
    }
        return fCtrs;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::swirlNumber::write()
{
    if (active_)
    {
        // NOTE: this is the main part of the function object and implements the
        // actual functionality.

        // Retrieve a reference to the velocity field
        const volVectorField& U = obr_.lookupObject<volVectorField>(velocityFieldName_);


	// Modified by Rahul: 18/10/21
        const volScalarField& rho =
            obr_.lookupObject<volScalarField>("rho");

	 const volScalarField& pd =
              obr_.lookupObject<volScalarField>("pd");
	// End of modification


        // interpolate onto the faces
	surfaceVectorField UFace = fvc::interpolate(U);
	surfaceScalarField rhoFace = fvc::interpolate(rho); // Rahul: 19/10/21
	surfaceScalarField pdFace = fvc::interpolate(pd);   // Rahul: 19/10/21

        // scalar distanceVal;

        //- Unit vectors to construct mutually orthogonal unit vectors to patchNormalVector
        // used in calculating tangential velocity component as OpenFOAM doesn't have builtin
        // implementation for cylindrical coordiante system. Calculations may seem trivial for 
        // faceZones whose normals are aligned with coordinate vectors.
        vector tangentUnitVector1(vector::zero);
        vector tangentUnitVector2(vector::zero);
        vector axialUnitVector(vector::zero);

        vectorField meshPoints = mesh_.Cf();

        vectorField zoneFaceCenters = faceCentres(faces_, meshPoints); //returning coordinate of the face center of facezone of mesh_
        
        for(label pointI = 0; pointI != zoneFaceCenters.size(); pointI++)
        {
            // distanceVal = mag(zoneFaceCenters[pointI] - originVector_);
            
            //- Mutually orthogonal unit vectors in faceZone
            if(mag(zoneFaceCenters[pointI] - originVector_) >= patchRadius_/5.) 
            {
                axialUnitVector = patchNormalVector_/mag(patchNormalVector_);
                
                tangentUnitVector1 = zoneFaceCenters[pointI] - originVector_; //radial in cylinderical coordinate
                tangentUnitVector1 /= mag(tangentUnitVector1);

                tangentUnitVector2 = (tangentUnitVector1 ^ axialUnitVector); //^ cross, perpendicular to  (radial and axial) - tangential in cylinderical coordinate
                tangentUnitVector2 /= mag(tangentUnitVector2);

                break; // break out of the loop as soon as unit vectors are determined 
            }
        }

        //- Calculation of integrals necessary for calculation of swirl number
        scalar GPhi = 0;
        scalar Gx = 0;
	

        if(!(faces_.size()==zoneFaceCenters.size()))
        {
            Info << "Problem running utility" << endl;
        }
        else
        {
            // for(label indxI = 0; indxI!=faces_.size(); indxI++)
    	   
            surfaceScalarField magTangentialU
                (mag((UFace&tangentUnitVector1)*tangentUnitVector1 + (UFace&tangentUnitVector2)*tangentUnitVector2));     
	 
            surfaceScalarField magAxialU(UFace&axialUnitVector);

            scalarField radiusField(mag(zoneFaceCenters - originVector_));
   

            forAll(faces_, indxI)
            {
		  // Modified by Rahul: 18/10/21
   	             GPhi += rhoFace[faces_[indxI]]*sqr(radiusField[indxI])*magAxialU[faces_[indxI]]*magTangentialU[faces_[indxI]]*mesh_.magSf()[faces_[indxI]];		 
	             Gx += (rhoFace[faces_[indxI]]*sqr(magAxialU[faces_[indxI]])+(pdFace[faces_[indxI]]-p_inf))*radiusField[indxI]*mesh_.magSf()[faces_[indxI]]; // pd at mixingTube-combustor center (at 0.03s for case 3_... - Syngas) = 110437                   
		  // End of modification
            }
            const scalar v_small = 1e-10;
        
        reduce(GPhi, sumOp<scalar>()); 
        reduce(Gx, sumOp<scalar>()); 
        
        
        if (Pstream::master())
       {
           // Call the base class method which checks if the output file exists
           // and creates it, if necessary. That also calls the .writeFileHeader()
           // method of the derived class.

           const scalar swirlNumber_ = GPhi/(patchRadius_*Gx + v_small); 
           
           logFiles::write();  

           Info << " Swirl number = " << swirlNumber_ << nl;
           
           Info << " Numerator number = " << GPhi << nl;
           Info << " Denomenator number = " << Gx << nl;
           // Add the entry for this time step that has just been computed.
           file(MAIN_FILE) << obr_.time().value() << tab << swirlNumber_ << endl;
           // writeTime(file(MAIN_FILE));
           // file(MAIN_FILE) << tab << swirlNumber_ << endl;
       }

        }
    } // 
    return true;
}
// ************************************************************************* //
