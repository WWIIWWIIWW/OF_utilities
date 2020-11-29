const bool multiComponent(thermo.lookupOrDefault<bool>("multiComponent", false));
const bool nonUnitySc(thermo.lookupOrDefault<bool>("nonUnitySc", false));

scalar Sct  = thermo.lookupOrDefault<scalar>("Sct", 1.0);
dictionary ScDict = thermo.subDict("Sc");
dictionary DijDict = thermo.subDict("Dij");

wordList speciesNames = composition.species();

List<scalar> Sc(Y.size());
PtrList<volScalarField> D_(Y.size());
PtrList<volScalarField> x_(Y.size());

forAll(speciesNames, i)
{
   Sc[i] = readScalar(ScDict.lookup(speciesNames[i]));
}

forAll(Y, i)
{
    x_.set
    (
        i, new volScalarField
        (
            IOobject
            (
                "x_" + speciesNames[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless
        )
    );

    D_.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "D_" + speciesNames[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("D_", dimMass / dimLength / dimTime, 0.0)
        )
    );
}

List<scalar> column(Y.size());
List<List<scalar> > Dij_(Y.size(), column);

forAll(Dij_,i)
{
    List<scalar>& tDij_ = Dij_[i];
    List<scalar> s(DijDict.lookup(speciesNames[i]));
    tDij_ = s;
}

volScalarField zero
(
    IOobject
    (
        "zero",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,-2,1,0,0,0,0), 0.0)
); 
//***********************************************************//
const scalar R = constant::physicoChemical::R.value()*1000;
if (multiComponent && !nonUnitySc)
{
    Info<< "MultiComponent Library ON!" << endl;
}
if (!multiComponent && nonUnitySc)
{
    Info<< "nonUnitySc ON!" << endl;
}
if (multiComponent && nonUnitySc)
{
    Info<< "MultiComponent Library ON and nonUnitySc ON is now allowed!" << endl;
    Info<< "Exit!" << endl;
    return 0;
}
else
{
    Info<< "Running in original reactingLMKsgsFoam" << nl << endl;
}
