{
    Info<< "Pre-calculating E & F" << endl;

    dimensionedScalar small("small", dimVelocity, SMALL);
    volScalarField uD = mag(U);
    volScalarField tuD = mag(fvc::laplacian(1e-6*(fvc::curl(U))));

    uD.ref() = tuD * mesh.V() *2.0e6;
    volScalarField s = (uD+small)/SL;

    Ka = Foam::pow(s,1.5);
    Ka.ref() = Ka * Foam::pow(dL*dL*dL/mesh.V(),1./6.);

    F.ref() = max(1.0, 5.* Foam::pow(mesh.V(),1./3.) / dL); 

    volScalarField alph = 2.0/(1.0+s)+2.35/(1.0+1/s)-2.0;
    E.ref() = max(1.0, Foam::pow(s,alph));
}
