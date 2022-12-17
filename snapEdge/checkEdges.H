
//label nOutline = meshOutline.size();

forAll(meshOutline, i)
{

    label fi0 = edgeBelongsToFeature[i].d0.featurei;
    label fi1 = edgeBelongsToFeature[i].d1.featurei;
    
    if ( (fi0 != -1) || (fi1 != -1) )
    {

	label i0 = edgeBelongsToFeature[i].d0.pointi;
	scalar di0 = edgeBelongsToFeature[i].d0.distance;
	vector pi0 = edgeBelongsToFeature[i].d0.p;
	label i1 = edgeBelongsToFeature[i].d1.pointi;
	scalar di1 = edgeBelongsToFeature[i].d1.distance;
	vector pi1 = edgeBelongsToFeature[i].d1.p;

	if (fi0 == -1)
	{
	    i0 = edgeBelongsToFeature[i].d1.pointo;
	    pi0 = newPoints[i0];
	    //di0 = edgeBelongsToFeature[i].d1.distanceo;
	    di0 = distance(globalStlPoints[is], stlFeatures[is][fi1], pi0);
	}

	if (fi1 == -1)
	{
	    i1 = edgeBelongsToFeature[i].d0.pointo;
	    pi1 = newPoints[i1];
	    //di1 = edgeBelongsToFeature[i].d0.distanceo;
	    di1 = distance(globalStlPoints[is], stlFeatures[is][fi0], pi1);
	}

	forAll(meshOutline, j)
	{
  	    if ( i != j )
	    {
		label fj0 = edgeBelongsToFeature[j].d0.featurei;
		label fj1 = edgeBelongsToFeature[j].d1.featurei;

		if ( ((fj0 != -1) || (fj1 != -1)) && ( (edgeBelongsToFeature[i].d0.featurei != -1) || (edgeBelongsToFeature[i].d1.featurei != -1)) )
		{
		  
		    label j0 = edgeBelongsToFeature[j].d0.pointi;
		    vector pj0 = edgeBelongsToFeature[j].d0.p;
		    scalar dj0 = edgeBelongsToFeature[j].d0.distance;
		    
		    label j1 = edgeBelongsToFeature[j].d1.pointi;
		    vector pj1 = edgeBelongsToFeature[j].d1.p;
		    scalar dj1 = edgeBelongsToFeature[j].d1.distance;

		    if (fj0 == -1)
		    {
			j0 = edgeBelongsToFeature[j].d1.pointo;
			pj0 = newPoints[j0];
			dj0 = edgeBelongsToFeature[j].d1.distanceo;
			dj0 = distance(globalStlPoints[is], stlFeatures[is][fj1], pj0);
		    }

		    if (fj1 == -1)
		    {
			j1 = edgeBelongsToFeature[j].d0.pointo;
			pj1 = newPoints[j1];
			dj1 = edgeBelongsToFeature[j].d0.distanceo;
			dj1 = distance(globalStlPoints[is], stlFeatures[is][fj0], pj1);
		    }

		    bool connected = (i0 == j0) || (i0 == j1) || (i1 == j0) || (i1 == j1);

		    if ( !connected )
		    {

			scalar lambda1, lambda2;
			vector px = vector::zero;
			//scalar scale = mag(pj1-pj0);
			scalar d1 = distance(pj0, pj1, pi0, lambda1, px);
			scalar d2 = distance(pj0, pj1, pi1, lambda2, px);

			bool over1 = (lambda1 >= minFit) && (lambda1 <= maxFit) && (d1 < overlap);
			bool over2 = (lambda2 >= minFit) && (lambda2 <= maxFit) && (d2 < overlap);

			if (over1 || over2)
			{
			    if (max(di0,di1) > max(dj0,dj1))
			    {
				edgeBelongsToFeature[i].d0.featurei = -1;
				edgeBelongsToFeature[i].d1.featurei = -1;
			    }
			    else
			    {
				edgeBelongsToFeature[j].d0.featurei = -1;
				edgeBelongsToFeature[j].d1.featurei = -1;
			    }
			} // if (overlapping)
		    }
		    else
		    {
		        // if the edges are connected, check if the points that are not connected are mapped onto each other
		      // disconnect the point on the line with the greatest angle
		        vector pi01 = pi1 - pi0;
			vector pj01 = pj1 - pj0;
			if ((j1 == i0) || (j0 == i1) )
			{
			    pj01 *= -1.0;
			}
			scalar dot = pi01 & pj01;
			if (dot > SMALL)
			{
			    // i0 is connected
			    if ((i0 == j0) || (i0 == j1))
			    {
			        scalar ai = edgeBelongsToFeature[i].d0.cos;
				if (i0 == j0)
				{
				    scalar aj = edgeBelongsToFeature[j].d0.cos;
				    if ( (ai>cosParAngle) && (aj>cosParAngle) )
				    {
					if (dj1 > di1)
					{
					    edgeBelongsToFeature[j].d1.featurei = -1;
					}
					else
					{
					    edgeBelongsToFeature[i].d1.featurei = -1;
					}
				    }
				    else
				    {
					if (ai > aj)
					{
					    edgeBelongsToFeature[j].d1.featurei = -1;
					}
					else
					{
					    edgeBelongsToFeature[i].d1.featurei = -1;
					}
				    }
				}
				else
				{
				    scalar aj = edgeBelongsToFeature[j].d1.cos;
				    if ( (ai>cosParAngle) && (aj>cosParAngle) )
				    {
					if (dj0 > di1)
					{
					    edgeBelongsToFeature[j].d0.featurei = -1;
					}
					else
					{
					    edgeBelongsToFeature[i].d1.featurei = -1;
					}
				    }
				    else
				    {
					if (ai > aj)
					{
					    edgeBelongsToFeature[j].d0.featurei = -1;
					}
					else
					{
					    edgeBelongsToFeature[i].d1.featurei = -1;
					}
				    }
				}
			    }

			    // i1 is connected
			    if ((i1 == j0) || (i1 == j1))
			    {
			        scalar ai = edgeBelongsToFeature[i].d1.cos;
				if (i1 == j0)
				{
				    scalar aj = edgeBelongsToFeature[j].d0.cos;
				    if ( (ai>cosParAngle) && (aj>cosParAngle) )
				    {
					if (dj1 > di0)
					{
					    edgeBelongsToFeature[j].d1.featurei = -1;
					}
					else
					{
					    edgeBelongsToFeature[i].d0.featurei = -1;
					}
				    }
				    else
				    {
					if (ai > aj)
					{
					    edgeBelongsToFeature[j].d1.featurei = -1;
					}
					else
					{
					    edgeBelongsToFeature[i].d0.featurei = -1;
					}
				    }
				}
				else
				{
				    scalar aj = edgeBelongsToFeature[j].d1.cos;
				    if ( (ai>cosParAngle) && (aj>cosParAngle) )
				    {
					if (dj0 > di0)
					{
					    edgeBelongsToFeature[j].d0.featurei = -1;
					}
					else
					{
					    edgeBelongsToFeature[i].d0.featurei = -1;
					}
				    }
				    else
				    {
					if (ai > aj)
					{
					    edgeBelongsToFeature[j].d0.featurei = -1;
					}
					else
					{
					    edgeBelongsToFeature[i].d0.featurei = -1;
					}
				    }
				}
			    }

			}
		    } // if (!connected)
		}
	    } // end if ( (i != j) && ...
	}
    }

}

