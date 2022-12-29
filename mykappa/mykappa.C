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

\*---------------------------------------------------------------------------*/

#include "mykappa.H"
#include "turbulenceModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"


#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "psiReactionThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "wallDist.H"
#include "boundBox.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mykappa, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        mykappa,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::mykappa::writeFileHeader(const label i)
{
    writeHeader(file(), "mykappa ()");

    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::mykappa::calcmykappa
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& tau_chemRR,
    const volScalarField& delta_
)
{
    //volScalarField y = wallDist(mesh_).y();
    //volScalarField center = mesh_.C().component(vector::X);
    
    //volScalarField y_ = y - center;
    
    tmp<volScalarField> tmykappa
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("tmykappa", dimless, 0.0)
        )
    );

    
    volScalarField& mykappa = tmykappa.ref();

    //tmp<volScalarField> ttau_mix(magSqr(y_) * rho / alpha);
    tmp<volScalarField> ttau_mix(magSqr(delta_) * rho / alpha);
    const volScalarField& tau_mix = ttau_mix();
    
    mykappa = 1 / (1+ (tau_mix / (tau_chemRR) + ROOTVSMALL));
 

    return tmykappa;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mykappa::mykappa
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log)
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::mykappa::~mykappa()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mykappa::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    return true;
}


bool Foam::functionObjects::mykappa::execute()
{
    if (mesh_.foundObject<compressible::turbulenceModel>(turbulenceModel::propertiesName))
    {
        const compressible::turbulenceModel& turbModel = mesh_.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );
        
        const psiReactionThermo& thermo = mesh_.lookupObject<psiReactionThermo>
        (
            psiReactionThermo::dictName
        );
        
        const volScalarField& tau_chemRR =
                     mesh_.lookupObject<volScalarField>("tau_chemRR");
                     
        #include "calcCellSize.H"
        
        word name(type());

        return store
        (
            name, 
            calcmykappa(turbModel.alphaEff(), thermo.rho(), tau_chemRR, delta_)
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::mykappa::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    return true;
}


// ************************************************************************* //
