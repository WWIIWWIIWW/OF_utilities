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
    reaction->correct();
    Qdot = reaction->Qdot();
    volScalarField Yt(0.0*Y[0]);

    if (multiComponent && !nonUnitySc) //multiComponent
    {
        #include "massDiffusive.H"
        forAll(Y, i)
        {
            if (i != inertIndex && composition.active(i))
            {
                volScalarField& Yi = Y[i];

                fvScalarMatrix YiEqn
                (
                    fvm::ddt(rho, Yi)
                  + mvConvection->fvmDiv(phi, Yi)
                  - fvm::laplacian(turbulence->mut()/Sct + D_[i], Yi)
                 ==
                    reaction->R(Yi)
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
    }
    if (!multiComponent && nonUnitySc) //nonUnitySc
    {
        forAll(Y, i)
        {
            if (i != inertIndex && composition.active(i))
            {
                volScalarField& Yi = Y[i];

                fvScalarMatrix YiEqn
                (
                    fvm::ddt(rho, Yi)
                  + mvConvection->fvmDiv(phi, Yi)
                  - fvm::laplacian(turbulence->mut()/Sct +turbulence->mu()/Sc[i], Yi)
                 ==
                    reaction->R(Yi)
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
    }

    else //Original
    {
        forAll(Y, i)
        {
            if (i != inertIndex && composition.active(i))
            {
                volScalarField& Yi = Y[i];

                fvScalarMatrix YiEqn
                (
                    fvm::ddt(rho, Yi)
                  + mvConvection->fvmDiv(phi, Yi)
                  - fvm::laplacian(turbulence->muEff(), Yi)
                 ==
                    reaction->R(Yi)
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
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
