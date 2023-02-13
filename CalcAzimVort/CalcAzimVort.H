volScalarField AzimVort
(
    IOobject
    (
        "AzimVort",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("AzimVort", dimVelocity / dimLength, 0.0)
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
    	
volVectorField radialUnitVector
(
    IOobject
    (
        "radialUnitVector",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("radialUnitVector", dimless, vector(1, 1, 1))
);
