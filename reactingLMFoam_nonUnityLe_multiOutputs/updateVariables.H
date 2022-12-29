//update variables; note in creatFieldREFs, we have only claimed variables
forAll(Y, i)
{
    ROP_[i].field() = -(reaction->R(Y[i])()).source()/mesh.V(); //ROP.Yi
    rho_[i].field() = rho * Y[i]; //rho.Yi
}

ksgs_       = turbulence->k(); //Instantaneous subgrid-scale tke
epsilonsgs_ = turbulence->epsilon();
musgs_      = turbulence->mu();
nusgs_      = turbulence->nu();
alphasgs_   = turbulence->alpha();
Cp_         = thermo.Cp();
tau_chemRR  = chem.tc();
#include "LESResIndex.H"

//write below if not AUTO_WRITE
//if (runTime.outputTime())
//{
//    ksgs_.write();
//    epsilonsgs_.write();
//    musgs_.write();
//    nusgs_.write();
//    alphasgs_.write();
//    Cp_.write();
//    //LESResIndex.write();
//}
//else
//{
//    Info << "Extra variables specified not written" << endl;
//}




