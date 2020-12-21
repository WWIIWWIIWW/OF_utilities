const dictionary& subChemProp = chemistryProperties.subDict("TFMcoeff");

dimensionedScalar SL("SL", dimVelocity, readScalar(subChemProp.lookup("SL")));
dimensionedScalar dL("dL", dimLength, readScalar(subChemProp.lookup("dL")));
dimensionedScalar maxQdot("maxQdot", dimEnergy/dimTime/dimVolume, readScalar(subChemProp.lookup("maxQdot")));

