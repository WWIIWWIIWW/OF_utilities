dictionary fuel(mixtureFractionDict.subDict("fuel"));
dictionary oxidizer(mixtureFractionDict.subDict("oxidizer"));
const word fractionBasis(mixtureFractionDict.lookup("fractionBasis"));

speciesTable species;
autoPtr<chemistryReader<gasHThermoPhysics>> pchemistry
(
    chemistryReader<gasHThermoPhysics>::New(thermo, species)
);
chemistryReader<gasHThermoPhysics>& chemistry = pchemistry();
const speciesCompositionTable& specieComposition = chemistry.specieComposition();

if ((fractionBasis != "mass") && (fractionBasis != "mole"))
{
    FatalError << "Unknown fractionBasis type"
               << token::SPACE
               << fractionBasis
               << nl
               << "Valid types are: mass or mole."
               << abort(FatalError);
}

scalarField c_Frac(Y.size());
scalarField h_Frac(Y.size());
scalarField o_Frac(Y.size());

scalarField Y_fuel(Y.size(), 0.0);
scalarField Y_oxi(Y.size(), 0.0);
scalarField X_fuel(Y.size(), 0.0);
scalarField X_oxi(Y.size(), 0.0);

scalar Yc_fuel, Yc_oxi;
scalar Yh_fuel, Yh_oxi;
scalar Yo_fuel, Yo_oxi;

Yc_fuel = 0;
Yc_oxi = 0;
Yh_fuel = 0;
Yh_oxi = 0;
Yo_fuel = 0;
Yo_oxi = 0;

forAll(species, specieI)
{
    scalar C = 0;
    scalar H = 0;
    scalar O = 0;
      
    forAll(specieComposition[species[specieI]], elemI)
    {
        word name = specieComposition[species[specieI]][elemI].name();
        label n = specieComposition[species[specieI]][elemI].nAtoms();
            
        if (name == "C")
        {
            C = atomicWeights["C"]*n/composition.W(specieI);
        }
        if (name == "H")
        {
            H = atomicWeights["H"]*n/composition.W(specieI);
        }
        if (name == "O")
        {
            O = atomicWeights["O"]*n/composition.W(specieI);
        }
    }

    c_Frac[specieI] = C;
    h_Frac[specieI] = H;
    o_Frac[specieI] = O;
}

scalar mwf = 0.0;
scalar mwo = 0.0;
if (fractionBasis == "mole")
{
    const scalar mTotf = sum(X_fuel);
    const scalar mToto = sum(X_oxi);
        
    forAll(species, specieI)
    {
        if (fuel.found(species[specieI]))
        {
            X_fuel[specieI] = readScalar(fuel.lookup(species[specieI]));
        }
        if (oxidizer.found(species[specieI]))
        {
            X_oxi[specieI] = readScalar(oxidizer.lookup(species[specieI]));
        } 

        X_fuel[specieI] /= mTotf;
        X_oxi[specieI] /= mToto;
            
        mwf += composition.W(specieI)*X_fuel[specieI];
        mwo += composition.W(specieI)*X_oxi[specieI];
    }
        
    forAll(species, specieI)
    {
        Y_fuel[specieI] = X_fuel[specieI]*composition.W(specieI)/mwf;
        Y_oxi[specieI] = X_oxi[specieI]*composition.W(specieI)/mwo;
    }
}
else
{
    forAll(species, specieI)
    {
        if (fuel.found(species[specieI]))
        {
            Y_fuel[specieI] = readScalar(fuel.lookup(species[specieI]));
        }
        if (oxidizer.found(species[specieI]))
        {
            Y_oxi[specieI] = readScalar(oxidizer.lookup(species[specieI]));
        }
    }
}

forAll(species, specieI)
{    
    Yc_fuel += Y_fuel[specieI]*c_Frac[specieI];
    Yh_fuel += Y_fuel[specieI]*h_Frac[specieI];
    Yo_fuel += Y_fuel[specieI]*o_Frac[specieI];

    Yc_oxi += Y_oxi[specieI]*c_Frac[specieI];
    Yh_oxi += Y_oxi[specieI]*h_Frac[specieI];
    Yo_oxi += Y_oxi[specieI]*o_Frac[specieI];
}
