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
    RemoveVortons

Description
    RemoveVortons saved in U -> due to the use of LEMOS by kai.zhang.1@city.ac.ul
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    OFstream file("time_list.txt");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        //file << runTime.timeName() << tab << endl;

        Foam::string src = runTime.timePath()/"U";
        Foam::string dst = runTime.timePath()/"U.bak";

        Foam::string cmd = "sed '/vortons/,/;/d' " + dst + ">" + src;

        if (fileHandler().isFile(dst.c_str()))
        {
            Info << "Vortons already removed for time "
                 << runTime.timeName() 
                 << ", Skiped!"
                 << endl;
        }
        else
        {
            Info<< "Removing vortorns from U, and create backup for time "
                << runTime.timeName() 
                << endl;

            cp(src.c_str(), dst.c_str());
            system(cmd.c_str());
        }
    }


    Info<< "\nEnd\n" << endl;
    
    return 0;
}


// ************************************************************************* //
