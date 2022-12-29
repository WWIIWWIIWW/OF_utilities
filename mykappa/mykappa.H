/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::functionObjects::mykappa

Description (under development, the version does not apply to other cases except the jet one we have)
    Evaluates and outputs kappa = (1/(1+Da)) 
    - Da = tau_mix/tau_c, if tau_mix<<tau_c, kappa~1, we well resolve near-wall chemistry.
    - tau_mix = mixing time scale = (cell size in wall normal direction)^2 / effective thermal diffusivity. (or mass diffusivty if Preff = 1)
    -tau_c = chemical time scale based on tc in:
        https://www.openfoam.com/documentation/guides/latest/api/PaSR_8C_source.html
    
    Values written to time directories as field 'mykappa'.

    Example of function object specification:
    \verbatim
    kappa
    {
        type            mykappa;
        libs            ( "libkappa.so" );
        ...
    }
    \endverbatim

Usage
    \table
        Property | Description                | Required   | Default value
        type     | type name: mykappa           | yes        |
    \endtable
    
    // wallDist dic must be added to fvSchemes (not necessary, 
    but the code was planned to calculate wallDist, so used))
    wallDist
    {
    method          meshWave;

        // Optionally correct distance from near-wall cells to the boundary
        correctWalls    true;
    }

Note
    Writing field 'mykappa' is done by default, but it can be overridden by
    defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    mykappa.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_mykappa_H
#define functionObjects_mykappa_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"

#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class mykappa Declaration
\*---------------------------------------------------------------------------*/

class mykappa
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate y+
        tmp<volScalarField> calcmykappa
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volScalarField& tau_chemRR,
            const volScalarField& delta_
        );

        //- Disallow default bitwise copy construct
        mykappa(const mykappa&);

        //- Disallow default bitwise assignment
        void operator=(const mykappa&);


public:

    //- Runtime type information
    TypeName("mykappa");


    // Constructors

        //- Construct from Time and dictionary
        mykappa
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~mykappa();


    // Member Functions

        //- Read the mykappa data
        virtual bool read(const dictionary&);

        //- Calculate the mykappa field
        virtual bool execute();

        //- Write the mykappa field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
