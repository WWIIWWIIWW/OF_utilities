#include "label.H"
#include "scalar.H"
#include "vector.H"

using namespace Foam;

#ifndef TYPES_H
#define TYPES_H

struct eDetails
{
    label pointi;
    label pointj;
    label feature;
    scalar distance;
    bool fixed;
    scalar cos;
    vector pi;
    vector pj;
};

struct details
{
    label featurei;
    label pointi;
    label pointo;
    scalar distance;
    scalar distanceo;
    scalar lambda;
    scalar lambdaOther;
    scalar cos;
    vector p;
    vector po;
};

struct edgeDetails
{
    details d0;
    details d1;
};

#endif
