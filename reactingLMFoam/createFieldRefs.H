const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
const label inertIndex(composition.species()[inertSpecie]);

Info<< "Create fields for ROP_.\n" << endl;
PtrList<volScalarField> ROP_(Y.size());
forAll(Y, i)
{
    const word RRname = "ROP." + Y[i].name();
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
}

