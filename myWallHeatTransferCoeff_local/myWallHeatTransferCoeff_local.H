/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::functionObjects::myWallHeatTransferCoeff_local

Description
    Calculates and write the HTC based on a fixed reference temperature

    \f[
            h = \frac{q}{T_{ref} - T_W}
    \f]

    at wall patches as the volScalarField field 'myWallHeatTransferCoeff_local'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    myWallHeatTransferCoeff_local1
    {
        type        myWallHeatTransferCoeff_local;
        libs        ("libfieldFunctionObjects.so");
        ...
        region      fluid;
        patches     (".*Wall");
    }
    \endverbatim

Usage
    \table
        Property | Description                        | Required   | Default value
        type     | type name: myWallHeatTransferCoeff_local | yes     
        patches  | list of patches to process         | no         | all wall patches
        region   | region to be evaluated             | no         | default region
    \endtable

Note
    Writing field 'myWallHeatTransferCoeff_local' is done by default, but it can be overridden by
    defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    myWallHeatTransferCoeff_local.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_myWallHeatTransferCoeff_local_H
#define functionObjects_myWallHeatTransferCoeff_local_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "HashSet.H"
#include "volFieldsFwd.H"

#include "wallFvPatch.H" //new

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class myWallHeatTransferCoeff_local Declaration
\*---------------------------------------------------------------------------*/

class myWallHeatTransferCoeff_local
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private Data

    //- Optional Tref (default is 0)
    //scalar Tref;

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the HTC
        tmp<volScalarField> calcmyWallHeatTransferCoeff_local
        (
            const volScalarField& alpha,
            const volScalarField& he,
            const volScalarField& T
        );


public:

    //- Runtime type information
    TypeName("myWallHeatTransferCoeff_local");


    // Constructors

        //- Construct from Time and dictionary
        myWallHeatTransferCoeff_local
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        myWallHeatTransferCoeff_local(const myWallHeatTransferCoeff_local&) = delete;


    //- Destructor
    virtual ~myWallHeatTransferCoeff_local();


    // Member Functions

        //- Read the myWallHeatTransferCoeff_local data
        virtual bool read(const dictionary&);

        //- Calculate the wall HTC
        virtual bool execute();

        //- Write the wall HTC
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const myWallHeatTransferCoeff_local&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
