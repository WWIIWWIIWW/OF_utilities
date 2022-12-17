#ifndef FUNCTIONS_H
#define FUNCTIONS_H

template<class T>
void addToList(List<T>& el, const T& e)
{
    label n = el.size();
    el.setSize(n+1);
    el[n] = e;
}

template<class T>
void addToList(List<T>& el, const List<T>& e)
{
    label n = el.size();
    label m = e.size();
    el.setSize(n+m);
    for(label i=0; i<m; i++)
    {
        el[n+i] = e[i];
    }
}

// calculates the distance for point p to line a-b 
// and returns the projection of p on line a-b in px
// lambda is the normalized distance of px from point a 
scalar distance
(
    const vector& a,
    const vector& b,
    const vector& p,
    scalar& lambda,
    vector& px
)
{
  /*
  Info << "a = " << a << endl;
  Info << "b = " << b << endl;
  Info << "p = " << p << ", mag(v) = " << mag(b-a) << endl;
  */
    vector v = b - a;
    scalar d = mag(p-a);
    lambda = 0.0;
    px = a;
    if (mag(v) > SMALL)
    {
        //Info << "hey stupido, its the same point" << endl;
        vector X = p - a;
	lambda = (X & v)/(v & v);
	px = a + lambda*v;
	d = mag(p - px);
    }
    return d;
}

// calculates the distance for point p to line a-b 
scalar distance
(
    const vector& a,
    const vector& b,
    const vector& p
)
{
    scalar lam = GREAT;
    vector px = vector::zero;
    return distance(a, b, p, lam, px);
}

// calculate the distance from a point to a line
scalar distance
(
    const vectorField& ps,
    const edge& e,
    const vector& p,
    scalar& lambda,
    vector& px
)
{
    vector a = ps[e[0]];
    vector b = ps[e[1]];
  
    return distance(a, b, p, lambda, px);
}

// calculate the distance from a point to a line
scalar distance
(
    const List<vector>& ps,
    const edge& e,
    const vector& p
)
{
    vector a = ps[e[0]];
    vector b = ps[e[1]];
  
    return distance(a, b, p);
}


bool insideTriangle
(
    const vector& p0,
    const vector& p1,
    const vector& p2,
    const vector& x
)
{
    vector Sf = (p1 - p0) ^ (p2 - p0);
    Sf /=mag(Sf);
    
    vector n0 = (p1 - p0) ^ Sf;
    vector n1 = (p2 - p1) ^ Sf;
    vector n2 = (p0 - p2) ^ Sf;

    scalar a0 = (x - p0) & n0;
    scalar a1 = (x - p1) & n1;
    scalar a2 = (x - p2) & n2;
    
    return (a0 <= 0) && (a1 <= 0) && (a2 <= 0);
}

/*
bool checkOverlap
(
    const List<vector>& points,
    const List<bool>& isEdgePoint,
    const vector& p,
    const label n,
    const scalar dist
)
{
    label np=points.size();
    bool overlap = false;
    label i=0;
    while ( (!overlap) && (i<np) )
    {
        if (isEdgePoint[i])
        {
            if ( (mag(points[i]-p) < dist) && (i != n) )
            {
                overlap = true;
            }
        }
        i++;
    }

    return overlap;
}
*/

// idea!!!
// create a triangle of the two connected edged.
// we can calculate useful properties

#endif
