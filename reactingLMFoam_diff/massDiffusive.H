Info << "Update Molar Fractions!" << endl;
forAll(Y, i)
{
    const scalarField& psiCells = psi.internalField();
    const scalarField& TCells = T.internalField();
    const scalarField& yCells = Y[i].internalField();

    scalarField& xCells = x_[i].ref();

    forAll(x_[i].internalField(), cellI)
    {
       xCells[cellI] = yCells[cellI] * psiCells[cellI] * TCells[cellI] * R / composition.W(i);
    }

    forAll(x_[i].boundaryField(), patchI)
    {
        const fvPatchScalarField& ppsi = psi.boundaryField()[patchI];
        const fvPatchScalarField& pT = T.boundaryField()[patchI];
        const fvPatchScalarField& py = Y[i].boundaryField()[patchI];

        fvPatchScalarField& px = x_[i].boundaryFieldRef()[patchI];

        forAll(px, faceI)
        {
            px[faceI] = py[faceI] * ppsi[faceI] * pT[faceI] * R / composition.W(i);
        }
    }
}

Info << "Update Binary Diffusion!" << endl;
forAll(Y, i)
{
    if (i != inertIndex && composition.active(i))
    {
        volScalarField tmpGamma = 0 * zero ;

        forAll(Y, j)
        {
            if (j != i)
            {
                tmpGamma += x_[j] / dimensionedScalar("Dij", dimensionSet(0,2,-1,0,0), Dij_[i][j]);
            }
        }

        D_[i] = rho * (1-x_[i]) / (tmpGamma + dimensionedScalar("SMALL", dimensionSet(0,-2,1,0,0), Foam::SMALL));

    }
}
