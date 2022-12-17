if (U.db().foundObject<volVectorField>("UMean") * U.db().foundObject<volScalarField>("pMean") * U.db().foundObject<volScalarField>("TMean") == 1) 
{
    //*******************************Preparation*********************************
    ksgs_  = turbulence->k(); //Instantaneous subgrid-scale tke m2/s2
    epsgs_ = turbulence->epsilon(); //Instantaneous subgrid-scale dissipation rate m2/s3
    mu_    = turbulence->mu(); //dynamic viscosity nu * rho, also used in Yi equation kg/(m.s)
    nu_    = turbulence->nu(); //kinematic viscosity m2/s
    rho    = thermo.rho(); //density kg/m3
    Cp_    = thermo.Cp(); //Specific heat capacity
    alpha  = turbulence->alpha();

    Pt     = mu_/alpha
    //Le     = preset 
    //Sc     = 
    //Nu

    const objectRegistry& db = U.db();

    tmp<volVectorField> tUPrime(U - db.lookupObject<volVectorField>("UMean"));
    const volVectorField& UPrime = tUPrime(); //U'

    tmp<volScalarField> tpPrime(p - db.lookupObject<volScalarField>("pMean"));
    const volScalarField& pPrime = tpPrime(); //p'

    tmp<volScalarField> tTPrime(T - db.lookupObject<volScalarField>("TMean"));
    const volScalarField& TPrime = tTPrime(); //T'

    tmp<volTensorField> tUPrimeUPrime(UPrime * UPrime);
    const volTensorField& UPrimeUPrime = tUPrimeUPrime(); //U'U'

    tmp<volVectorField> tUPrimeTPrime(UPrime * TPrime);
    //const volVectorField& UPrimeTPrime = tUPrimeTPrime(); //U'T'
    UPrimeTPrime = tUPrimeTPrime();

    tmp<volTensorField> tgradUPrime(fvc::grad(UPrime));
    const volTensorField& gradUPrime = tgradUPrime(); //grad(U')

    tmp<volSymmTensorField> tSPrime(symm(gradUPrime));
    const volSymmTensorField& SPrime = tSPrime(); //symm(grad(U'))

    tmp<volVectorField> tgradpPrime(fvc::grad(pPrime));
    const volVectorField& gradpPrime = tgradpPrime(); //grad(p')

    tmp<volVectorField> tgradTPrime(fvc::grad(TPrime));
    const volVectorField& gradTPrime = tgradTPrime (); //grad(T')

    //tmp<volSymmTensorField> tTSPrime(twoSymm(gradUPrime));
    //const volSymmTensorField& TSPrime = tTSPrime(); //twosymm(grad(U'))
 
    tmp<volTensorField> tsigma(2 * mu * SPrime - 2.0/3.0 * mu * fvc::div(UPrime));
    const volTensorField& sigma = tsigma(); //sigma, applied to dissipation of kinetic energy, if i=!j, div(UPrime) = 0
    //***************************************************************************



    if (runTime.outputTime())
    {
        //patch1_U.write();
        //patch2_U.write();
        //#include "DNSbudget.H"
        //#include "write.H"
        //#include "write_T.H"
        #include "writeBudgets.H"
    }
}

else
{
    Info << "Running reactingDNS without budget!" << endl;
}
