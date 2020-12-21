rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, pd)); // UNSURE IF pd OR p HERE!!!

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(pd, rho, U, phiHbyA, rhorAUf, MRF); // UNSURE IF pd OR p HERE!!!

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pdEqn
    (
        fvc::ddt(psi) * pReff
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, pd)
     ==
        fvOptions(psi, pd, rho.name())
    );

    pdEqn.solve(mesh.solver(pd.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + pdEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
pd.relax();

U = HbyA - rAU*fvc::grad(pd);
U.correctBoundaryConditions();
fvOptions.correct(U);

/*
if (pressureControl.limit(pd))
{
    pd.correctBoundaryConditions();
    rho = thermo.rho();
}
*/