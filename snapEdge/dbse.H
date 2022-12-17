/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

Class
    Foam::dbse

Description

SourceFiles
    dbseI.H
    dbse.C
    dbseIO.C

\*---------------------------------------------------------------------------*/

#ifndef dbse_H
#define dbse_H

#include "fvMesh.H"
#include "types.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
//class Istream;
//class Ostream;

// Forward declaration of friend functions and operators
//class dbse;
  //Istream& operator>>(Istream&, dbse&);
  //Ostream& operator<<(Ostream&, const dbse&);


/*---------------------------------------------------------------------------*\
                          Class dbse Declaration
\*---------------------------------------------------------------------------*/

class dbse
{

    // Private data
        Time& runTime_;
        fvMesh& mesh_;
        dictionary dict_;

    // some useful references
        const polyBoundaryMesh& patches_;
        const faceZoneMesh& zones_;
        const pointField& points_;

  
        wordList snapPatches_;
        wordList snapZones_;
        wordList stlFileNames_;

        scalar tolerance_;
        scalar relaxation_;

        scalar featureAngle_;
        scalar excludeAngle_;
        scalar parallelAngle_;

        bool includeInterior_;
        label nIterations_;

        scalar fitFactor_;
        scalar minFit_;
        scalar maxFit_;

        vectorField newPoints_;
        labelListList edgePoints_;
        label nSTL_;
        List<edgeList> stlFeatures_;
        List<vectorField> globalSTLPoints_;
        List< List<label> > nEdgePoints_;

        scalar cosFeature_;
        scalar cosFeature2_;
        scalar cosExclude_;
        scalar cosParAngle_;

        scalar smallestEdgeLength_;

public:

    // Static data members
        const static scalar deg2rad_;

    // Constructors
        //- Construct from components
        dbse
        (
	    Time& runTime,
	    fvMesh& mesh,
	    const dictionary& snapEdgeDict
        );

    //- Destructor
        ~dbse();


    // Member Functions

        void snapEdges();

        void traceMethod
        (
	    const edgeList& meshOutline,
	    const labelList& addr,
	    const labelListList& pointEdges,
	    const labelList& edgeOutlineConnection
	);

        label findBestEdge(const List<eDetails>& edgeConnections);

        void traverseEdges
	(
	    const vectorField& points,
	    const edgeList& stlFeatures,
	    List<eDetails>& edgeConnections,
	    const List<List<label> >& edgeEdges,
	    const label& edgei,
	    label& nCounter
	);

        void connectEdgeToFeature
	(
	    const vectorField& points,
	    const edgeList& stlFeatures,
	    const label& i0,
	    const label& i1,
	    eDetails& edgeConnection
	);

        void recalculate
	(
	    const vectorField& points,
	    const edgeList& stlFeatures,
	    const label& edgei,
	    List<eDetails>& edgeConnections
	);

        List<label> findConnectedEdges
	(
	    const edgeList& meshOutline,
	    const labelList& addr,
	    const labelListList& pointEdges,
	    const labelList& edgeOutlineConnection,
	    const label& edgei
	);

        bool edgesSharePoint
        (
	    const edge& e1,
	    const edge& e2
	);

        const vectorField& newPoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "dbseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
