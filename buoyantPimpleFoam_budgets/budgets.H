if (U.db().foundObject<volScalarField>("rhoMean") * U.db().foundObject<volScalarField>("pMean") * U.db().foundObject<volVectorField>("rhoUMean") * U.db().foundObject<volScalarField>("rhoTMean") == 1) 
{
    //*******************************Preparation*********************************
    const objectRegistry& db = U.db();

    //tmp<volScalarField> tpMean(db.lookupObject<volScalarField>("pMean"));
    //const volScalarField& pMean = tpMean(); //rhoMean

    tmp<volScalarField> trhoMean(db.lookupObject<volScalarField>("rhoMean"));
    const volScalarField& rhoMean = trhoMean(); //rhoMean

    tmp<volVectorField> trhoUMean(db.lookupObject<volVectorField>("rhoUMean"));
    const volVectorField& rhoUMean = trhoUMean(); //rhoUMean

    tmp<volScalarField> trhoTMean(db.lookupObject<volScalarField>("rhoTMean"));
    const volScalarField& rhoTMean = trhoTMean(); //rhoTMean

    //*******************************favre**************************************
    FavUMean = rhoUMean/rhoMean;
    FavTMean = rhoTMean/rhoMean;

    FavUPrime = U - rhoUMean/rhoMean;
    FavTPrime = T - rhoTMean/rhoMean;

    FavUPrimeUPrime = FavUPrime * FavUPrime;
    FavUPrimeTPrime = FavUPrime * FavTPrime;
    FavTPrimeTPrime = FavTPrime * FavTPrime;

    FavrhoUPrimeUPrime = rho * FavUPrimeUPrime;
    FavrhoUPrimeTPrime = rho * FavUPrimeTPrime;
    FavrhoTPrimeTPrime = rho * FavTPrimeTPrime;

    tmp<volTensorField> tFavgradUPrime(fvc::grad(FavUPrime));
    const volTensorField& FavgradUPrime = tFavgradUPrime();

    tmp<volSymmTensorField> tFavTSPrime(twoSymm(FavgradUPrime));
    const volSymmTensorField& FavTSPrime = tFavTSPrime();

    const tensor I(1, 0, 0, 0, 1, 0, 0, 0, 1);
    tmp<volTensorField> tFavsigma(nu_ * (FavTSPrime - 2.0/3.0 * I * tr(FavgradUPrime)));
    const volTensorField& Favsigma = tFavsigma ();

    epsilon_ = Favsigma && FavgradUPrime;
    k_ = 0.5 * magSqr(FavUPrime);

    //Integral time scale
    tau_mixI = k_/(epsilon_+epsilonSmall);
    //Kolmogorov time scale
    tau_mixK = Foam::sqrt(mag(nu_/(epsilon_+epsilonSmall)));
    //Geometric mean of Kolmogorov and integral time scales
    tau_mixM = Foam::sqrt(mag(tau_mixI * tau_mixK));

    //Integral length scale
    L_I = Foam::pow(k_, 3.0 / 2.0) / epsilon_;
    //Kolmogorov length scale
    //volScalarField a = ;
    L_K = Foam::pow(Foam::pow(nu_, 3.0) / epsilon_, 1.0/4.0);
    //kappa calculations
    R_IK = L_I / (L_K + lengthSmall); //ratio between integral length scale and kolmogorov scale
    R_IM = R_IK;
    R_IM.ref() = L_I / Foam::pow(mesh.V(), 1.0/3.0); //ratio between integral length scale and average mesh size

    if (runTime.outputTime())
    {
        #include "writeBudgets.H"
    }
    else
    {
        Info << "runTime.outputTime():" << runTime.outputTime() << endl;
    }
}

else
{
    Info << "Running reactingDNS without budget!" << endl;
}

    //*******************************Post-pro**************************************
    /*
    if (U.db().foundObject<volVectorField>("FavrhoTPrimeTPrimeMean") * U.db().foundObject<volScalarField>("FavrhoUPrimeTPrimeMean") * U.db().foundObject<volScalarField>("FavrhoTPrimeTPrimeMean") == 1)
    {
        Info << "S, B terms not written.28" << endl;
        tmp<volTensorField> tFavrhoUPrimeUPrimeMean (db.lookupObject<volTensorField>("FavrhoTPrimeTPrimeMean"));
        const volTensorField& FavrhoUPrimeUPrimeMean = tFavrhoUPrimeUPrimeMean (); //read mean of favre averaged rhoU'U'
        Info << "S, B terms not written.29" << endl;
        tmp<volVectorField> tFavrhoUPrimeTPrimeMean (db.lookupObject<volVectorField>("FavrhoUPrimeTPrimeMean"));
        const volVectorField& FavrhoUPrimeTPrimeMean = tFavrhoUPrimeTPrimeMean (); //read mean of favre averaged rhoU'T'
        Info << "S, B terms not written.30" << endl;
        tmp<volScalarField> tFavrhoTPrimeTPrimeMean (db.lookupObject<volScalarField>("FavrhoTPrimeTPrimeMean"));
        const volScalarField& FavrhoTPrimeTPrimeMean = tFavrhoTPrimeTPrimeMean (); //read mean of favre averaged rhoT'T'
        Info << "S, B terms not written.31" << endl;
        tmp<volScalarField> tS(FavrhoUPrimeUPrimeMean && fvc::grad(FavUMean));
        const volScalarField S = tS; //shear contribution to turbuent kinetic energy
        Info << "S, B terms not written.32" << endl;
        tmp<volScalarField> tB(FavUMean && fvc::grad(pMean));
        const volScalarField& B = tB (); //buoyancy contribution to turbuent kinetic energy
        Info << "S, B terms not written.33" << endl;
        if (runTime.outputTime())
        {
            S.write();
            B.write();
        }
        else
        {
            Info << "S, B terms not written.34" << endl;
        }
    }

    else
    {
        Info << "Running reactingDNS without post-processing."  << endl;
    }
    */

    //***************************************************************************
    /*stress
    tensor I(1, 0, 0, 0, 1, 0, 0, 0, 1);
    tmp<volTensorField> tFavsigma(nu_ * (FavTSPrime - 2.0/3.0 * I * tr(FavgradUPrime)));
    const volTensorField& Favsigma = tFavsigma ();

    //Epsilon, kinetic energy dissipation rate
    epsilon_      = FavgradUPrime.component(tensor::XX) * Favsigma.component(tensor::XX)+FavgradUPrime.component(tensor::YY) * Favsigma.component(tensor::YY) + FavgradUPrime.component(tensor::ZZ) * Favsigma.component(tensor::ZZ);
    */


    /*
    forAll(epsilon_, i)
    {
        if (epsilon_[i] > 0)
        {
            kappaIFR[i] = tau_chemFR[i] / (tau_chemFR[i] + tau_mixI[i] + timeSmall[i]);
            kappaKFR[i] = tau_chemFR[i] / (tau_chemFR[i] + tau_mixK[i] + timeSmall[i]);
            kappaMFR[i] = tau_chemFR[i] / (tau_chemFR[i] + tau_mixM[i] + timeSmall[i]);
            kappaIRR[i] = tau_chemRR[i] / (tau_chemRR[i] + tau_mixI[i] + timeSmall[i]);
            kappaKRR[i] = tau_chemRR[i] / (tau_chemRR[i] + tau_mixK[i] + timeSmall[i]);
            kappaMRR[i] = tau_chemRR[i] / (tau_chemRR[i] + tau_mixM[i] + timeSmall[i]);
        }
        else
        {
            kappaIFR[i] = 1.0;
            kappaKFR[i] = 1.0;
            kappaMFR[i] = 1.0;
            kappaIRR[i] = 1.0;
            kappaKRR[i] = 1.0;
            kappaMRR[i] = 1.0;
        }
    }
    */    



