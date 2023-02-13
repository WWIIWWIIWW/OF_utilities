/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ODE2

Description
    $ wmake
    $ ODE2 RKCK45 > log
    ODE solver for solving mass-spring system with intial displacement,
    with explanation/comments by Kai Zhang
Ref: http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2008/ZongyuanGu/reportZongyuan.pdf
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//ODEsystem is an abstract class defines the system of first order ODEs as explained above. 
//ODEsolver is the base class for all the ODE solvers in OpenFOAM.
class myODE:public ODESystem {
    private:
        //- Mass of the system
        const scalar m_;
        //- Stiffness of the system
        const scalar k_;
        
    public:
        //constructor function
        //let's call parent class constructor
        //and initialise the child class private variables
        myODE(const scalar& mass, const scalar& stiffness):ODESystem(),m_(mass),k_(stiffness){}

        //Member functions
        label nEqns() const
        {
            return 2;
        }

        void derivatives
        (
            const scalar t, //was x
            const scalarField& y,
            scalarField& dydt //was dydx
        ) const
        {
            dydt[0] = y[1];
            dydt[1] = -(k_/m_+VSMALL)*y[0];
        }

        void jacobian
        (
            const scalar t,
            const scalarField& y,
            scalarField& dfdt,
            scalarSquareMatrix& dfdy
        ) const
        {
        /*
            dfdy[0] = 0.0;
            dfdy[1] = 0.0;

            dfdy[0][0] = 0.0;
            dfdy[0][1] = 1.0;

            dfdy[1][0] = -(k_/m_+VSMALL);
            dfdy[1][1] = 0;
            */
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
//argc (argument count) and argv (argument vector)
int main(int argc, char *argv[])
{
    //allow the program to read the ODEsolver name from the command line. 
    argList::validArgs.append("ODESolver");
    argList args(argc, argv);
    
    // system properties (mass and stiffness)
    const scalar m = 1.0; // Kg
    const scalar k = 1.0; // N/m
    
    // Initial displacement
    const scalar y0 = 2e-03; // m
    const label N = 200;
    
    //time setting
    scalar dtEst = 0.01; // integration initial step , must be smaller than the actual time step ->check my chemstry setting in other solvers 
    scalar startTime = 0.0; // s
    const scalar endTime = 100; // s
    const scalar dt = endTime/N; //this is the time step for other solvers
    
    // Create the instance of myODE
    myODE ode(m,k);

    // Create dictionary and add the odeSolver name
    dictionary dict;
    dict.add("solver", args[1]);

    // Create the selected ODE system solver with (ode, solverName)
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict); //New is a seletor pointing to 

    // Initialise the ODE system fields
    // Initial displacement and velocity
    
    scalarField yStart(ode.nEqns());
    yStart[0] = y0;
    yStart[1] = 0.0;


    // Required to store dydx
    scalarField dydt(ode.nEqns()); //initialise dydt to be sizeof(nEqns), used to store constructed dydt

    // Integration loop
    scalar tEnd;
    for (label t=0; t<N; t++)
    {
        tEnd = startTime + dt; //integrate by time step
        ode.derivatives(startTime, yStart, dydt); //(t, y, dydt)
        odeSolver->solve(startTime, tEnd, yStart, dtEst);
        startTime = tEnd; //reset startTime to the end time point
        Info << startTime << "   " << yStart[0] << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
