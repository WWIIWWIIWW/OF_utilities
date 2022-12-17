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

\*---------------------------------------------------------------------------*/

#include "dbse.H"
#include "functions.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::dbse::staticData();
const Foam::scalar Foam::dbse::deg2rad_ = M_PI/180.0;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dbse::dbse
(
    Time& runTime,
    fvMesh& mesh,
    const dictionary& snapEdgeDict
)
:
    runTime_(runTime),
    mesh_(mesh),
    dict_(snapEdgeDict),
    patches_(mesh.boundaryMesh()),
    zones_(mesh.faceZones()),
    points_(mesh.points()),
    snapPatches_(dict_.lookup("snapPatches")),
    snapZones_(dict_.lookup("snapZones")),
    stlFileNames_(dict_.lookup("stlFileNames")),
    tolerance_(readScalar(dict_.lookup("tolerance"))),
    relaxation_(readScalar(dict_.lookup("relaxation"))),
    featureAngle_(readScalar(dict_.lookup("featureAngle"))),
    excludeAngle_(readScalar(dict_.lookup("excludeEdgeAngle"))),
    parallelAngle_(readScalar(dict_.lookup("parallelAngle"))),
    includeInterior_(readBool(dict_.lookup("includeInterior"))),
    nIterations_(readLabel(dict_.lookup("nIterations"))),
    fitFactor_(readScalar(dict_.lookup("fitFactor"))),
    minFit_(-mag(fitFactor_)),
    maxFit_(1.0 + mag(fitFactor_)),

    newPoints_(points_),
    edgePoints_(0),
    nSTL_(stlFileNames_.size()),
    stlFeatures_(nSTL_),
    globalSTLPoints_(nSTL_),
    nEdgePoints_(nSTL_),

    cosFeature_(::cos(featureAngle_*deg2rad_)),
    cosFeature2_(::cos(featureAngle_*0.5*deg2rad_)),
    cosExclude_(::cos(excludeAngle_*deg2rad_)),
    cosParAngle_(::cos(parallelAngle_*deg2rad_)),
    smallestEdgeLength_(GREAT)

{

    forAll(stlFileNames_, is)
    {
        word stlFileName(stlFileNames_[is]);
	fileName stlFile("constant/triSurface/"+stlFileName);
        triSurface stlSurface(stlFile);
        
        const edgeList& stlEdges = stlSurface.edges();    
        const labelListList& stlEdgeFaces = stlSurface.edgeFaces();
        const vectorField& stlSf = stlSurface.faceNormals();
        const vectorField& stlPoints = stlSurface.localPoints();
        addToList(globalSTLPoints_[is], stlPoints);
        nEdgePoints_[is].setSize(stlPoints.size());

        Info << "Finding features for stl : " << stlFileName << endl;
        forAll(stlEdgeFaces, i)
        {
            if(!stlSurface.isInternalEdge(i))
            {
                addToList(stlFeatures_[is], stlEdges[i]);
                nEdgePoints_[is][stlEdges[i][0]]++;
                nEdgePoints_[is][stlEdges[i][1]]++;
            }
            else
            {
                label fi0 = stlEdgeFaces[i][0];
                label fi1 = stlEdgeFaces[i][1];
        
                scalar magF0 = mag(stlSf[fi0]);
                scalar magF1 = mag(stlSf[fi1]);
                if ( (magF0 < 1.0e-5) || (magF1 < 1.0e-5))
                {
                    Info << "Warning:: face normals are zero!?!" << endl;
                }
                vector f0 = stlSf[fi0]/magF0;
                vector f1 = stlSf[fi1]/magF1;
        
                scalar cosa = f0 & f1;
                if (cosa < cosFeature_)
                {
                    addToList(stlFeatures_[is], stlEdges[i]);
                    nEdgePoints_[is][stlEdges[i][0]]++;
                    nEdgePoints_[is][stlEdges[i][1]]++;
                }
            }
        }

        if (stlFeatures_[is].size() == 0)
        {
            Info << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
            Info << "WARNING!!! Your stl does not contain any feature lines." << endl;
            Info << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
        }

    }

    Info << "Done!" << endl;

    // check that all patches are available and get the smallest length scale of the mesh
    forAll(snapPatches_, index)
    {

        label id = patches_.findPatchID(snapPatches_[index]);
        if (id == -1)
        {
            Info << "Could not find patch: " << snapPatches_[index] <<endl;
            Info << "Available patches are:" << endl;
            forAll(patches_, i)
            {
                Info << patches_[i].name() << endl;
            }
            FatalError << abort(FatalError);
        }
        else
        {
            const labelList& addr = patches_[id].meshPoints();
            const edgeList& edges = patches_[id].edges();
        
            forAll(edges, i)
            {
                const edge& e = edges[i];
                label i0 = addr[e[0]];
                label i1 = addr[e[1]];
                scalar scale = mag(points_[i0] - points_[i1]);
                smallestEdgeLength_ = min(smallestEdgeLength_, scale);
            }
        }
    }

    // simple copy/paste/adapt for faceZones
    forAll(snapZones_, index)
    {
        label id = zones_.findZoneID(snapZones_[index]);
        if (id == -1)
        {
            Info << "Could not find zone: " << snapZones_[index] <<endl;
            Info << "Available zones are:" << endl;
            Info << zones_.names() << endl;
            FatalError << abort(FatalError);
        }
        else
        {
            const labelList& addr = zones_[id]().meshPoints();
            const edgeList& edges = zones_[id]().edges();
        
            forAll(edges, i)
            {
                const edge& e = edges[i];
                label i0 = addr[e[0]];
                label i1 = addr[e[1]];
                scalar scale = mag(points_[i0] - points_[i1]);
                smallestEdgeLength_ = min(smallestEdgeLength_, scale);
            }
        }
    }

    Info << "smallestEdgeLength = " << smallestEdgeLength_ << endl;

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dbse::~dbse()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::dbse::snapEdges()
{

    for(label ni=0; ni<nIterations_; ni++)
    {
        vectorField oldNewPoints(newPoints_);

        Info << "(" << ni+1 << "/" << nIterations_ << ") Matching edges...";
        flush(Info);

        forAll(snapPatches_, index)
        {

            label id = patches_.findPatchID(snapPatches_[index]);

            const labelList& addr = patches_[id].meshPoints();
            const vectorField& bSf = patches_[id].faceNormals();
            const vectorField& bCf = patches_[id].faceCentres();

            const edgeList& edges = patches_[id].edges();
            const labelListList& ef = patches_[id].edgeFaces();
            const labelListList& pointEdges = patches_[id].pointEdges();

            labelList edgeOutlineConnection(ef.size(), -1);
            edgeList meshOutline(0);

            // construct the list of all the mesh edges to be used for matching
            forAll(ef, edgei)
            {
                if (!patches_[id].isInternalEdge(edgei))
                {
                    edgeOutlineConnection[edgei] = meshOutline.size();
                    addToList(meshOutline, edges[edgei]);
                }
                else
                {
                    if (includeInterior_)
                    {
                        label fi0 = ef[edgei][0];
                        label fi1 = ef[edgei][1];
              
                        scalar mag0 = mag(bSf[fi0]);
                        scalar mag1 = mag(bSf[fi1]);
                        if ( (mag0<SMALL) || (mag1<SMALL) )
                        {
                            Info << "Warning:: zero area faces in the mesh!!!" << endl;
                        }
                        vector f0 = bSf[fi0]/(mag0 + VSMALL);
                        vector f1 = bSf[fi1]/(mag1 + VSMALL);
                        vector v1 = bCf[fi1] - bCf[fi0];
                        v1 /= mag(v1) + VSMALL;
                        if ( (f0 & v1) > 0)
                        {
                            f0 *= -1.0;
                        }

                        if ( (f1 & v1) < 0)
                        {
                            f1 *= -1.0;
                        }
                        // check angle between faces (use a smaller feature angle than the one used for the stl)
                        scalar cosa = mag(f0 & f1);
                        if (cosa < cosFeature2_)
                        {
                            edgeOutlineConnection[edgei] = meshOutline.size();
                            addToList(meshOutline, edges[edgei]);
                        }
                    }
                }
            }

	    traceMethod(meshOutline, addr, pointEdges, edgeOutlineConnection);

        }

        // again: simple copy/paste/adapt for faceZones
        forAll(snapZones_, index)
        {

            label id = zones_.findZoneID(snapZones_[index]);

            const labelList& addr = zones_[id]().meshPoints();
            const vectorField& bSf = zones_[id]().faceNormals();
            const edgeList& edges = zones_[id]().edges();
            const labelListList& ef = zones_[id]().edgeFaces();
            const labelListList& pointEdges = patches_[id].pointEdges();
            labelList edgeOutlineConnection(edges.size(), -1);
            edgeList meshOutline(0);
        
            // construct the list of all the mesh edges to be used for matching
            forAll(ef, i)
            {
                if (!zones_[id]().isInternalEdge(i))
                {
                    edgeOutlineConnection[i] = meshOutline.size();
                    addToList(meshOutline, edges[i]);
                }
                else
                {
                    if (includeInterior_)
                    {
                        label fi0 = ef[i][0];
                        label fi1 = ef[i][1];
            
                        scalar mag0 = mag(bSf[fi0]);
                        scalar mag1 = mag(bSf[fi1]);
                        if ( (mag0<SMALL) || (mag1<SMALL) )
                        {
                            Info << "Warning:: zero area faces in the mesh!!!" << endl;
                        }
                        vector f0 = bSf[fi0]/(mag0 + VSMALL);
                        vector f1 = bSf[fi1]/(mag1 + VSMALL);
            
                        // check angle between faces (use a smaller feature angle than the one used for the stl)
                        scalar cosa = mag(f0 & f1);
                        if (cosa < cosFeature2_)
                        {
                            edgeOutlineConnection[i] = meshOutline.size();
                            addToList(meshOutline, edges[i]);
                        }
                    }
                }
            }

	    traceMethod(meshOutline, addr, pointEdges, edgeOutlineConnection);

        }

        Info << "Moving points... ";
        flush(Info);
	Info << sum(mag(newPoints_-oldNewPoints)) << endl;

        mesh_.movePoints(newPoints_);

        Info << "Done!" << endl;

    }

}

void Foam::dbse::traceMethod
(
    const edgeList& meshOutline,
    const labelList& addr,
    const labelListList& pointEdges,
    const labelList& edgeOutlineConnection
)
{

    eDetails dummy;
    dummy.feature = -1;
    dummy.pointi = -1;
    dummy.pointj = -1;
    dummy.fixed = false;
    dummy.distance = GREAT;
    dummy.pi = vector::zero;
    dummy.pj = vector::zero;
    
    label nMeshOutline = meshOutline.size();
    
    for(label is=0; is<nSTL_; is++)
    {

        List<eDetails> edgeConnections(nMeshOutline, dummy);
        List< List<label> > edgeEdges(nMeshOutline);
	label nSize = edgeEdges.size();
	Info << "find edgeEdges (" << nSize << ")";
	flush(Info);
	label counter = 0;
	forAll(edgeEdges, i)
        {
            counter++;
	    label a = (counter*10)/nSize;
	    if (a == 1)
	    {
                counter = 0;
		Info << ".";
		flush(Info);
	    }
	    List<label> ce = findConnectedEdges(meshOutline, addr, pointEdges, edgeOutlineConnection, i);
	    edgeEdges[i].setSize(ce.size());
	    forAll(ce, j)
	    {
                edgeEdges[i][j] = ce[j];
	    }
	}

	Info << "Connecting edges to features...";
	// connect all edges to a feature
	for(label i=0; i<nMeshOutline; i++)
        {
            const edge& edgei = meshOutline[i];
   
            connectEdgeToFeature
            (
                globalSTLPoints_[is], 
                stlFeatures_[is],
		addr[edgei[0]],
		addr[edgei[1]],
		edgeConnections[i]
	    );
	}

        Info << "findBestEdge...";
	// find the best fitting edge
	label nCurrentEdge = findBestEdge(edgeConnections);
	label nCounter = 1;
	Info << "traversing...";
	flush(Info);
	while (nCurrentEdge != -1)
	{
	    // traverse all connected segments
	    traverseEdges
	    (
	        globalSTLPoints_[is],
		stlFeatures_[is],
		edgeConnections, 
		edgeEdges, 
		nCurrentEdge,
		nCounter
	    );
	    nCurrentEdge = findBestEdge(edgeConnections);
	}
	Info << "Done" << endl;
	flush(Info);

	forAll(edgeConnections, i)
	{
            bool fixed = edgeConnections[i].fixed;
	    bool assoc = (edgeConnections[i].feature != -1);
	    if ((!fixed) && assoc)
	    {
		  //Info << ". This edge is left hanging..." << i << endl;
	    }
	    else
	    {
	        if (assoc)
	        {
		    scalar dist = edgeConnections[i].distance;
		    if (dist > 1.0e-10)
		    {
			//Info << "edgei = " << i << ", dist = " << dist << endl;
		    }
		}
	    }
	    //Info << endl;
	}
    }  // end for(label is=0...

}

Foam::label Foam::dbse::findBestEdge
(
 const Foam::List<eDetails>& edgeConnections 
)
{

    label iClose = -1;
    scalar minDist = GREAT;
    for (label i=0; i<edgeConnections.size(); i++)
    {
        bool fixed = edgeConnections[i].fixed;
        if (!fixed)
        {
            scalar dist = edgeConnections[i].distance;
            if (dist < minDist)
            {
                minDist = dist;
                iClose = i;
            }
        }
    }

    return iClose;

}

void Foam::dbse::traverseEdges
(
    const vectorField& points,
    const edgeList& stlFeatures,
    List<eDetails>& edgeConnections,
    const List<List<label> >& edgeEdges,
    const label& edgei,
    label& nCounter
)
{
    label i0 = edgeConnections[edgei].pointi; 
    label i1 = edgeConnections[edgei].pointj;
    //Info << "edgei = " << edgei << ", i0 = " << i0 << ", i1 = " << i1 << endl;
    const List<label>& iEdges = edgeEdges[edgei];
    label nEdges = iEdges.size();

    // check for overlap before moving points on all fixed connected edges
    bool overlap = false;
    
    vector pi = edgeConnections[edgei].pi;
    vector pj = edgeConnections[edgei].pj;
    vector pAv = 0.5*(pi+pj);
    if (mag(pi-pj) < 0.1*smallestEdgeLength_)
    {
      //Info << "warning...edge is collapsing" << endl;
        label nEdgeEdges = edgeEdges[edgei].size();
	bool moved = false;
	label i = 0;
	while ( (i<nEdgeEdges) && (!moved) )
	{
	    label ie = edgeEdges[edgei][i];
	    label xi = edgeConnections[ie].pointi;
	    label xj = edgeConnections[ie].pointj;
	    if ( (i0 == xi) || (i0 == xj) || (i1 == xi) || (i1 == xj) )
	    {

	        if (edgeConnections[edgei].feature == edgeConnections[ie].feature)
		{
		    moved = true;

		    if (i0 == xi)
		    {
			vector pa = 0.5*(pj + edgeConnections[ie].pj);
			pi = pa;
			edgeConnections[edgei].pi = pi;
		    }
		    else
		    {
		        if (i0 == xj)
			{
			    vector pa = 0.5*(pj + edgeConnections[ie].pi);
			    pi = pa;
			    edgeConnections[edgei].pi = pi;
			}
			else
			{
			    if (i1 == xi)
			    {
				vector pa = 0.5*(pi + edgeConnections[ie].pj);
				pj = pa;
				edgeConnections[edgei].pj = pj;
			    }
			    else
			    {
				vector pa = 0.5*(pi + edgeConnections[ie].pi);
				pj = pa;
				edgeConnections[edgei].pj = pj;
			    }
			}
		    }
		}
	    }
	    i++;
	}
    }

    //scalar distFeat = edgeConnections[edgei].distance;
    for(label ie=0; ie<edgeConnections.size(); ie++)
    {
        if (ie == edgei) ie++;

        bool fixed = edgeConnections[ie].fixed;
	bool assoc = edgeConnections[ie].feature != -1;
        if (fixed && assoc)
        {
            label ii = edgeConnections[ie].pointi;
            label ij = edgeConnections[ie].pointj;

            vector pti = edgeConnections[ie].pi;
            vector ptj = edgeConnections[ie].pj;
	    vector ptAv = 0.5*(pti+ptj);

            scalar lam = GREAT;
            vector px = vector::zero;
            //scalar distanceToFeature = edgeConnections[ie].distance;

            if ( (i0 != ii) && (i0 != ij) )
            {
                scalar dist = distance(pti, ptj, pi, lam, px);
                if ((lam >= minFit_) && (lam <= maxFit_) && (dist < 1.0e-6))
                {
                    overlap = true;
                }
            }
            if ( (i1 != ii) && (i1 != ij) )
            {
                scalar dist = distance(pti, ptj, pj, lam, px);
                if ((lam >= minFit_) && (lam <= maxFit_) && (dist < 1.0e-6)) 
                {
                    overlap = true;
                }
            }

            if ( (ii != i0) && (ii != i1) )
            {
                scalar dist = distance(pi, pj, pti, lam, px);
                if ((lam >= minFit_) && (lam <= maxFit_) && (dist < 1.0e-6))
                {
                    overlap = true;
                }
            }
            if ( (ij != i0) && (ij != i1) )
            {
                scalar dist = distance(pi, pj, ptj, lam, px);
                if ((lam >= minFit_) && (lam <= maxFit_) && (dist < 1.0e-6)) 
                {
                    overlap = true;
                }
            }


        }
    }

    if (!overlap)
    {
        newPoints_[i0] += relaxation_*(pi - newPoints_[i0]);
        newPoints_[i1] += relaxation_*(pj - newPoints_[i1]);
    }

    edgeConnections[edgei].fixed = true;
    /*
       mesh.movePoints(newPoints);
       runTime++;
       mesh.write();
    */

    recalculate
    (
       points,
       stlFeatures,
       edgei,
       edgeConnections
    );

    // chose the best edge to continue with
    label iBest = -1;
    scalar criteria = GREAT;
    for(label i=0; i<nEdges; i++)
    {
        label ie = iEdges[i];
	bool assoc = edgeConnections[ie].feature != -1;
	if (!edgeConnections[ie].fixed && assoc)
	{
	    if (edgeConnections[ie].distance < criteria)
	    {
                //criteria = edgeConnections[ie].cos;
	        criteria = edgeConnections[ie].distance;
		iBest = ie;
	    }
	}
    }

    if (iBest >= 0)
    {
        nCounter++;
		
                //Info << "traverseEdges" << endl;
                traverseEdges
                (
                   points,
                   stlFeatures,
                   edgeConnections, 
                   edgeEdges, 
                   iBest,
                   nCounter
               );
    }
}

void Foam::dbse::connectEdgeToFeature
(
    const vectorField& points,
    const edgeList& stlFeatures,
    const label& i0,
    const label& i1,
    eDetails& edgeConnection
)
{

    vector p0 = newPoints_[i0];
    vector p1 = newPoints_[i1];
    scalar scale = mag(p1 - p0);
    vector n0 = (p1 - p0)/scale;

    edgeConnection.pointi = i0;
    edgeConnection.pointj = i1;

    forAll(stlFeatures, j)
    {
        const edge& feat = stlFeatures[j];
        vector gp0 = points[feat[0]];
        vector gp1 = points[feat[1]];
        vector n1 = gp1 - gp0;
        n1 /= mag(n1);
        
        scalar ac = mag(n1 & n0);
        
        // only align features which aren't close to perpendicular
        if (ac > cosExclude_)
        {
            scalar lam = GREAT;
            vector pxi = vector::zero;
	    vector pAv = 0.5*(p0 + p1);
            scalar di = distance(gp0, gp1, pAv, lam, pxi);
        
            if ( (lam >= minFit_) && (lam <= maxFit_) && ( di < tolerance_*scale) )
            {
                scalar lami = GREAT;
		scalar lamj = GREAT;
		vector p0i = vector::zero;
		vector p0j = vector::zero;
		distance(gp0, gp1, p0, lami, p0i);
		distance(gp0, gp1, p1, lamj, p0j);

                if (edgeConnection.feature == -1)
                {
                    edgeConnection.feature = j;
                    edgeConnection.distance = di;
                    edgeConnection.cos = ac;
		    edgeConnection.pi = p0i;
		    edgeConnection.pj = p0j;
                }
                else
                {
                    scalar acp = edgeConnection.cos;
                    // if both features are close to parallel, use distance, otherwise use angle
                    if ( (acp > cosParAngle_) && (ac > cosParAngle_) )
                    {
                        if (di < edgeConnection.distance)
                        {
                            edgeConnection.feature = j;
                            edgeConnection.distance = di;
                            edgeConnection.cos = ac;
			    edgeConnection.pi = p0i;
			    edgeConnection.pj = p0j;
                        }
                    }
                    else if ( ac > acp )
                    {
                        edgeConnection.feature = j;
                        edgeConnection.distance = di;
                        edgeConnection.cos = ac;
			edgeConnection.pi = p0i;
			edgeConnection.pj = p0j;
                    }
                }
            }
        
        }
    } // forAll(stlFeatures...
}

void Foam::dbse::recalculate
(
    const vectorField& points,
    const edgeList& stlFeatures,
    const label& edgei,
    List<eDetails>& edgeConnections
)
{
    label n = edgeConnections.size();
    label pi = edgeConnections[edgei].pointi;
    label pj = edgeConnections[edgei].pointj;

    for(label i=0; i<n; i++)
    {
        label ii = edgeConnections[i].pointi;
        label ij = edgeConnections[i].pointj;
        bool modified = false;

        if (!edgeConnections[i].fixed)
        {
            if ((ii==pi) || (ii==pj))
            {
                modified = true;
            }
            if ((ij==pi) || (ij==pj))
            {
                modified = true;
            }
        }

        if (modified)
        {
            connectEdgeToFeature
            (
                points,
                stlFeatures,
                ii,
                ij,
                edgeConnections[i]
            );
        }
    }
}

Foam::List<Foam::label> Foam::dbse::findConnectedEdges
(
    const edgeList& meshOutline,
    const labelList& addr,
    const labelListList& pointEdges,
    const labelList& edgeOutlineConnection,
    const label& edgei
)
{
    List<label> conEdges(0);

    label pi = meshOutline[edgei][0];
    label pj = meshOutline[edgei][1];

    labelList il = pointEdges[pi];
    for(label id=0; id<il.size(); id++)
    {
        label i = edgeOutlineConnection[il[id]];
    
        if (edgei != i)
        {
            bool connected = edgesSharePoint(meshOutline[edgei], meshOutline[i]);
        
            if ( connected )
            {
                label nSize = conEdges.size();
                conEdges.setSize(nSize+1);
                conEdges[nSize] = i;
            }
        }
    }
    
    labelList ij = pointEdges[pj];
    
    for(label id=0; id<ij.size(); id++)
    {
        label i = edgeOutlineConnection[ij[id]];
        if (edgei != i)
        {
            bool connected = edgesSharePoint(meshOutline[edgei], meshOutline[i]);
         
            if ( connected )
            {
                label nSize = conEdges.size();
                conEdges.setSize(nSize+1);
                conEdges[nSize] = i;
            }
        }
    }

    return conEdges;
}

bool Foam::dbse::edgesSharePoint
(
    const edge& e1,
    const edge& e2
)
{
    label i1 = e1[0];
    label j1 = e1[1];
    label i2 = e2[0];
    label j2 = e2[1];

    bool share = false;
    if ( (i1 == i2) || (i1 == j2) )
    {
        if ( i1 != -1)
        {
            share = true;
        }
    }
    if ( (j1 == i2) || (j1 == j2) ) 
    {
        if (j1 != -1)
        {
            share = true;
        }
    }

    return share;
}

const Foam::vectorField& Foam::dbse::newPoints()
{
    return newPoints_;
}

// ************************************************************************* //
