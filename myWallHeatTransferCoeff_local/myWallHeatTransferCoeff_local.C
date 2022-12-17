/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "myWallHeatTransferCoeff_local.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(myWallHeatTransferCoeff_local, 0);
    addToRunTimeSelectionTable(functionObject, myWallHeatTransferCoeff_local, dictionary);
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::myWallHeatTransferCoeff_local::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Wall HTC_local");
    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "integral");
    file() << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::myWallHeatTransferCoeff_local::calcmyWallHeatTransferCoeff_local
(
    const volScalarField& alpha,
    const volScalarField& he,
    const volScalarField& T
)
{
    tmp<volScalarField> tmyWallHeatTransferCoeff_local
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
            dimensionedScalar
            (
                "0",
                //dimMass/pow3(dimTime)/(dimTemperature/dimLength),
            	dimMass/(pow3(dimTime)*dimTemperature),
                0
            )
        )
    );
    
    tmp<volScalarField> tTref
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
            dimensionedScalar
            (
                "0",
            	dimTemperature,
                0
            )
        )
    );

    volScalarField::Boundary& myWallHeatTransferCoeff_localBf =
        tmyWallHeatTransferCoeff_local.ref().boundaryFieldRef();
        
    volScalarField::Boundary& TrefBF =
        tTref.ref().boundaryFieldRef();
        
    const volScalarField::Boundary& heBf =
        he.boundaryField();

    const volScalarField::Boundary& alphaBf =
        alpha.boundaryField();

    const volScalarField::Boundary& TBf =
        T.boundaryField();

    //find our celllist for next wall cells
    
    const fvPatchList& patches = mesh_.boundary();
    //labelList wallList;
    //wallList.clear();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];
        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                //get faceCelli for cell above the patch.
                label faceCelli = curPatch.faceCells()[facei];
                //wallList.resize(wallList.size()+1);
                //wallList[wallList.size()-1]=faceCelli;
                TrefBF[patchi][facei] = T[faceCelli];
            }
        }
    }
    
    
    forAll(myWallHeatTransferCoeff_localBf, patchi)
    {
        if (!myWallHeatTransferCoeff_localBf[patchi].coupled())
        {
            myWallHeatTransferCoeff_localBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad()/(TBf[patchi]-TrefBF[patchi]+ROOTVSMALL);
        }
    }

    if (foundObject<volScalarField>("qr"))
    {
        const volScalarField& qr = lookupObject<volScalarField>("qr");

        const volScalarField::Boundary& radHeatFluxBf =
            qr.boundaryField();

        forAll(myWallHeatTransferCoeff_localBf, patchi)
        {
            if (!myWallHeatTransferCoeff_localBf[patchi].coupled())
            {
                myWallHeatTransferCoeff_localBf[patchi]= (alphaBf[patchi]*heBf[patchi].snGrad()-radHeatFluxBf[patchi])/(TBf[patchi]-TrefBF[patchi]+ROOTVSMALL);
            }
        }
    }

    return tmyWallHeatTransferCoeff_local;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::myWallHeatTransferCoeff_local::myWallHeatTransferCoeff_local
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log),
    patchSet_()
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::myWallHeatTransferCoeff_local::~myWallHeatTransferCoeff_local()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::myWallHeatTransferCoeff_local::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall HTC on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::myWallHeatTransferCoeff_local::execute()
{
    word name(type());

    if
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        return store
        (
            name,
            calcmyWallHeatTransferCoeff_local(turbModel.alphaEff(), turbModel.transport().he(), turbModel.transport().T())
        );
    }
    else if (foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            lookupObject<solidThermo>(solidThermo::dictName);

        return store(name, calcmyWallHeatTransferCoeff_local(thermo.alpha(), thermo.he(), thermo.T()));
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::myWallHeatTransferCoeff_local::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& myWallHeatTransferCoeff_local =
        obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = myWallHeatTransferCoeff_local.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            file()
                << mesh_.time().value()
                << tab << pp.name()
                << tab << minHfp
                << tab << maxHfp
                << tab << integralHfp
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
