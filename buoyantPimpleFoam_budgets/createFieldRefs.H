const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();

volVectorField UPrimeTPrime
(
    IOobject
    (
        "UPrimeTPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("UPrimeTPrime", dimVelocity*dimTemperature, vector(0, 0, 0))
);


volTensorField UU
(
    IOobject
    (
        "UU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U*U
);

volVectorField UT
(
    IOobject
    (
        "UT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U*T
);

volScalarField TT
(
    IOobject
    (
        "TT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T*T
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoT
(
    IOobject
    (
        "rhoT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*T
);
volVectorField FavUMean
(
    IOobject
    (
        "FavUMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);
volScalarField FavTMean
(
    IOobject
    (
        "FavTMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T
);
volVectorField FavUPrime
(
    IOobject
    (
        "FavUPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);
volScalarField FavTPrime
(
    IOobject
    (
        "FavTPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T
);

volTensorField FavUPrimeUPrime
(
    IOobject
    (
        "FavUPrimeUPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U*U
);

volVectorField FavUPrimeTPrime
(
    IOobject
    (
        "FavUPrimeTPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U*T
);

volScalarField FavTPrimeTPrime
(
    IOobject
    (
        "FavTPrimeTPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T*T
);


volVectorField FavrhoUPrimeTPrime
(
    IOobject
    (
        "FavrhoUPrimeTPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U*T
);

volTensorField FavrhoUPrimeUPrime
(
    IOobject
    (
        "FavrhoUPrimeUPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U*U
);

volScalarField FavrhoTPrimeTPrime
(
    IOobject
    (
        "FavrhoTPrimeTPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*T*T
);


/*
volTensorField rhoUU
(
    IOobject
    (
        "rhoUU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U*U
);

volVectorField rhoUT
(
    IOobject
    (
        "rhoUT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U*T
);

volScalarField rhoTT
(
    IOobject
    (
        "rhoTT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*T*T
);
*/
volScalarField tau_mixI
(
    IOobject
    (
        "tau_mixI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tau_mixI", dimTime, 0.0)
);

volScalarField tau_mixK
(
    IOobject
    (
        "tau_mixK",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tau_mixK", dimTime, 0.0)
);

volScalarField tau_mixM
(
    IOobject
    (
        "tau_mixM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tau_mixM", dimTime, 0.0)
);

volScalarField L_I
(
    IOobject
    (
        "L_I",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("L_I", dimLength, 0.0)
);

volScalarField L_K
(
    IOobject
    (
        "L_K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("L_K", dimLength, 0.0)
);

volScalarField k_
(
    IOobject
    (
        "k_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    turbulence->k()
);

volScalarField epsilon_
(
    IOobject
    (
        "epsilon_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    turbulence->epsilon()
);

volScalarField mu_
(
    IOobject
    (
        "mu_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    turbulence->mu()
);

volScalarField nu_
(
    IOobject
    (
        "nu_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    turbulence->nu()
);

volScalarField Cp_
(
    IOobject
    (
        "Cp_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.Cp()
);

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    turbulence->alpha()
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.kappa()
);

volScalarField kappaIFR
(
    IOobject
    (
        "kappaIFR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaIFR", dimless, 0.0)
);

volScalarField kappaKFR
(
    IOobject
    (
        "kappaKFR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaKFR", dimless, 0.0)
);

volScalarField kappaMFR
(
    IOobject
    (
        "kappaMFR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaMFR", dimless, 0.0)
);

volScalarField kappaIRR
(
    IOobject
    (
        "kappaIRR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaIRR", dimless, 0.0)
);

volScalarField kappaKRR
(
    IOobject
    (
        "kappaKRR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaKRR", dimless, 0.0)
);

volScalarField kappaMRR
(
    IOobject
    (
        "kappaMRR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kappaMRR", dimless, 0.0)
);

const volScalarField epsilonSmall
(
    IOobject
    (
        "epsilonSmall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilonSmall", dimVelocity*dimVelocity/dimTime, SMALL)
);

volScalarField R_IK
(
    IOobject
    (
        "R_IK",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("R_IK", dimless, 0.0)
);

volScalarField R_IM
(
    IOobject
    (
        "R_IM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("R_IM", dimless, 0.0)
);

const volScalarField timeSmall
(
    IOobject
    (
        "timeSmall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("timeSmall", dimTime, SMALL)
);

const volScalarField lengthSmall
(
    IOobject
    (
        "lengthSmall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("lengthSmall", dimLength, SMALL)
);


