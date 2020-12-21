{
    Info<< "Pre-calculating E & F" << endl;
    //-----------------------------------------------//
    //----------Introducing the TFLES model----------//
    //-----------------------------------------------//
    /*
The summary version is that you need to replace instances of "boundaryField()" and "internalField()" with "boundaryFieldRef()" and "internalFieldRef()" if you want to modify them in your code. The non-"Ref" calls are const now.

Update

Seems like internalField and internalFieldRef have been replaced by direct calls to the field variable. Eg: U.internalField() becomes U() for a const reference to the internal field, or U.ref() for non-const access. See associated bug report.
         tmp<volScalarField> tepsilon(this->turbulence().epsilon());
         const scalarField& epsilon = tepsilon();
*/

    dimensionedScalar small("small", dimVelocity, SMALL);
    //O. Colin, F. Ducros, D. Veynante, T. Poinsot Phys. Fluids, 12 (7) (2000), pp. 1843-1863
    volScalarField uD = mag(U); //initialization;
    volScalarField tuD = mag(fvc::laplacian(1e-6*(fvc::curl(U))));
    uD.ref() = tuD * mesh.V() *2.0e6;

    //scalarField for following calculation
    volScalarField s = (uD+small)/SL;

    //Karlovitz number Ka
    Ka = Foam::pow(s,1.5);
    Ka.ref() = Ka * Foam::pow(dL*dL*dL/mesh.V(),1./6.);

    //Thickening factor F
    F.ref() = max(1.0, 5.* Foam::pow(mesh.V(),1./3.) / dL); 

    //Ratio kappa = E/F where E is the TFLES efficiency correction 
    volScalarField alph = 2.0/(1.0+s)+2.35/(1.0+1/s)-2.0;
    E.ref() = max(1.0, Foam::pow(s,alph));

}
