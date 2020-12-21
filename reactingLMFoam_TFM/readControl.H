const dictionary& subChemProp = chemistryProperties.subDict("TFMcoeff");
//Flame speed for the TFLES Model
dimensionedScalar SL("SL", dimVelocity, readScalar(subChemProp.lookup("SL")));
//Flame characteristic size for the TFLES Model
dimensionedScalar dL("dL", dimLength, readScalar(subChemProp.lookup("dL")));
//Maximum of heat release for the TFLES Model
dimensionedScalar maxQdot("maxQdot", dimEnergy/dimTime/dimVolume, readScalar(subChemProp.lookup("maxQdot")));

