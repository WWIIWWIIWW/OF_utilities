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

Class
    Foam::swirlNumber

Description
    Calculates and writes output of swirlNumber over a prescribed faceZone.

    \heading functionObject usage on a faceZone

    \table
        Property        | Description     | Required    | Default value
        faceZoneName    | name of facezone| yes         | -
        pointZoneName   | pointzone name  | yes         | -
        radius          | radius of zone  | yes         | 
        origin          | centerpoint     | yes         | (0 0 0)
        normal          | normal to zone  | yes         | (0 0 0)
    \endtable

    Example of the functionObject specification:
    \verbatim
    swirlNumberCalc
    {
        functionObjectLibs ("libswirlNumber.so");
        type swirlNumber;
        faceZoneName        planeFaceZone;
        pointZoneName       planePointZone;
        velocityFieldName   U;
        radius              0.02625;
        origin              (0 0 0);
        normal              (1 0 0);
    
        writeControl       timeStep;
        writeInterval      100;
    }
    \endverbatim

SourceFiles
    swirlNumber.H

\*---------------------------------------------------------------------------*/

#ifndef swirlNumber_H
#define swirlNumber_H

#include "volFieldsFwd.H"
#include "Switch.H"
#include "fvc.H"
#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "PstreamReduceOps.H"
#include "writeLocalObjects.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class swirlNumber Declaration
\*---------------------------------------------------------------------------*/

class swirlNumber
:
    public fvMeshFunctionObject,
    public logFiles
{
    protected:
    // Private data

        //- Name of the functionObject
        word name_;

        //- switch
        bool active_;

        //- Name of velocity field
        word velocityFieldName_;

        //- Name of the facezone
        word faceZoneName_;

        //- Name of the pointZone
        // word pointZoneName_;

        //- label of the facezone
        label faceZoneLabel_;

        //- label of the pointZone. Works best for pointZone created out of faceSet.
        // label pointZoneLabel_;
        
        //- radius of the marked faceZone
        scalar patchRadius_;

	// Rahul: 8/2/22: Adding p_inf
        //- Pressure of the fluid at rest (Ref: Litvinov 2021 - Swirl number and nozzle confinement effects in a flat-vane axial swirler)
        scalar p_inf;
        // End of modification
	
        //- normal vector to plane patch
        vector patchNormalVector_;

        //- origin vector on the patch any vector for reference will suffice
        // vector must lie on the patch
        vector originVector_;

        //- List of all labels
        // const labelList& faces_;
        const labelList& faces_; 
	
        //- labelList of all points
        // const labelList& points_;
        // const labelList& points_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        swirlNumber(const swirlNumber&);

        //- Disallow default bitwise assignment
        void operator=(const swirlNumber&);

    //protected:
         // NOTE: define a list of files this object writes to; the underlying logFiles
         // object will handle output to correct streams automatically.
         enum fileID
         {
             MAIN_FILE = 0
         };
     
         // NOTE: Create file names for each of the output files
         wordList createFileNames(const dictionary& dict) const;
     
         // NOTE: this first gets declared in the baseline object and gets used to
         // write the header in the output file.
         virtual void writeFileHeader(const label i);

public:


    //- Type information at runtime
    TypeName("swirlNumber");

    // Constructors

        //- Construct from given time and dictionary
        swirlNumber
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict 
        );

    // Destructor
    virtual ~swirlNumber();
    
    // Member Functions
 
        //- Return name of the functionObject in the controlDict
        virtual const word& name() const { return name_; }
    
        //- Read the swirlNumber data
        virtual bool read(const dictionary&);
    
        //- Execute, currently does nothing
        virtual bool execute();
    
        //- Execute at the final time-loop, currently does nothing
        virtual bool end();
    
        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();
    
        // NOTE: Key method which implements the object's actual functionality
        virtual bool write();
    
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&) {}
    
        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&) {}
    
    private:

        //- Calculation of face centers
        Foam::vectorField faceCentres(const labelList&, vectorField&);

        // Foam::surfaceVectorField faceCentres(const labelList&, vectorField&);
        
        Foam::surfaceScalarField calcTangentialU
            (surfaceVectorField&, labelList&, vector&, vector&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
