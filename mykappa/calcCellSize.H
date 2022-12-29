const faceList & ff = mesh_.faces();
const pointField & pp = mesh_.points();

volScalarField deltay = wallDist(mesh_).y();

forAll(mesh_.C(), celli)
{
    const cell& cc = mesh_.cells()[celli];
    labelList pLabels(cc.labels(ff));

    vectorField pLocal(pLabels.size(), vector::zero); 

    forAll(pLabels, pointi)
    {
        pLocal[pointi] = pp[pLabels[pointi]];
    }
    //dimensionedScalar Convert ("Convert", dimensionSet(0, 1, 0, 0, 0, 0, 0), scalar(1));
    deltay[celli] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));// * Convert;
}

volScalarField& delta_ = deltay;
delta_.ref() = mesh_.V() / magSqr(deltay.internalField());



  

