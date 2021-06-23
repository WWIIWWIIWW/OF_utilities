Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo
(
    psiReactionThermo::New(mesh)
);

psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

volScalarField YC
(
    IOobject
    (
        "YC",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
volScalarField YH
(
    IOobject
    (
        "YH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
volScalarField YO
(
    IOobject
    (
        "YO",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField Z
(
    IOobject
    (
        "Z",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Z", dimless, 0.0)
);
