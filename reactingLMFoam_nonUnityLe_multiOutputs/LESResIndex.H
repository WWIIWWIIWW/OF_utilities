//from https://github.com/WWIIWWIIWW/pimpleTKEBudgetFoam/blob/master/tkeBudget.H

if (U.db().foundObject<volVectorField>("UMean"))
{

    const objectRegistry& db = U.db();

    tmp<volVectorField> tUPrime(U - db.lookupObject<volVectorField>("UMean"));
    const volVectorField& UPrime = tUPrime(); //U'
    
    kRes =0.5*(UPrime & UPrime);//Instantaneous resolved tke
    kTot = kRes + ksgs_;
}

else
{
     Warning << "LESindex not calculated due to missing 'UMean' in database." << endl;
}
