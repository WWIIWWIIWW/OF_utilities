const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
const label inertIndex(composition.species()[inertSpecie]);

Info<< "Create fields for ROP_ and rho_.\n" << endl;
PtrList<volScalarField> ROP_(Y.size());
PtrList<volScalarField> rho_(Y.size());
forAll(Y, i)
{
    const word RRname  = "ROP." + Y[i].name();
    const word rhoname = "rho." + Y[i].name();
    
    ROP_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                RRname,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    ); 
    rho_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                rhoname,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimVolume, 0.0)
        )
    );  
}
//turbulence
volScalarField ksgs_
(
    IOobject
    (
        "ksgs_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->k()
);

volScalarField epsilonsgs_
(
    IOobject
    (
        "epsilonsgs_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    turbulence->epsilon()
);

volScalarField musgs_
(
    IOobject
    (
        "musgs_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->mu()
);

volScalarField nusgs_
(
    IOobject
    (
        "nusgs_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->nu()
);



volScalarField alphasgs_
(
    IOobject
    (
        "alphasgs_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->alpha()
);

//thermo
volScalarField Cp_
(
    IOobject
    (
        "Cp_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.Cp()
);

// LESindex
volScalarField kTot
(
    IOobject
    (
        "kTot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "kTot",
        dimVelocity*dimVelocity,
        0.0
    )
);

volScalarField kRes
(
    IOobject
    (
        "kRes",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "kRes",
        dimVelocity*dimVelocity,
        0.0
    )
);

//- ADDED CODE

volVectorField heatFluxCorr
(
    IOobject
    (
        "heatFluxCorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("heatFluxCorr", thermo.he().dimensions()/dimLength, Foam::vector(0,0,0))
);

volVectorField gradY
(
    IOobject
    (
        "gradY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradY", dimless/dimLength, Foam::vector(0,0,0))
);


dimensionedScalar Hs("Hs",thermo.he().dimensions(),0.0);
label Ysize = Y.size();

scalarList LewisNumber(Ysize,1.0);

Info<< "\nReading Lewis numbers\n" << endl;

dictionary LewisNumberDict(chemistryProperties.subDict("LewisNumbers"));


for (label i=0; i<Ysize; i++)
{
    if (LewisNumberDict.found(Y[i].member()))
    {
        LewisNumber[i] = readScalar(LewisNumberDict.lookup(Y[i].member()));
    }
    Info<<"Lewis number of specie "<<Y[i].name()<<" is: "<<LewisNumber[i]<<endl;
}


autoPtr<BasicChemistryModel<psiReactionThermo>> tempChem
(
    BasicChemistryModel<psiReactionThermo>::New(thermo)
);
BasicChemistryModel<psiReactionThermo>& chem = tempChem();

volScalarField tau_chemRR
(
    IOobject
    (
        "tau_chemRR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    chem.tc()
);

