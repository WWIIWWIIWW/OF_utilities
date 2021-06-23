Info<< "Calculating mixture fraction for time "
    << runTime.timeName() 
    << endl;

Yc.ref() *= 0;
Yh.ref() *= 0;
Yo.ref() *= 0;

forAll(species, specieI)
{
    Yc.ref() += Y[specieI]*c_Frac[specieI];
    Yh.ref() += Y[specieI]*h_Frac[specieI];
    Yo.ref() += Y[specieI]*o_Frac[specieI];
}

Z =
(
    (
        2.0*(Yc - Yc_oxi)/atomicWeights["C"]
      + (Yh - Yh_oxi)/(2.0*atomicWeights["H"])
      - (Yo - Yo_oxi)/atomicWeights["O"]
    )
   /(
        2.0*(Yc_fuel - Yc_oxi)/atomicWeights["C"]
      + (Yh_fuel - Yh_oxi)/(2.0*atomicWeights["H"])
      - (Yo_fuel - Yo_oxi)/atomicWeights["O"]
    )
);
