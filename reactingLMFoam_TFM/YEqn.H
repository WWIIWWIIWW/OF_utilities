tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    #include "calThickingFactor.H"
    reaction->correct();
    Qdot = reaction->Qdot();
    flame_sensor = tanh(100*mag(Qdot)/maxQdot);

    F = 1 + (F-1)*flame_sensor;
    E = 1 + (E-1)*flame_sensor;

    kappaTFM = E / F;

    Info<< "min/max(kappaTFM) = "
        << min(kappaTFM).value() << ", " << max(kappaTFM).value() << endl;

    Qdot *= E/F;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian((1-flame_sensor)*turbulence->muEff(), Yi)
              - fvm::laplacian(E*F*turbulence->alpha(), Yi) 
             ==
                E / F * reaction->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
