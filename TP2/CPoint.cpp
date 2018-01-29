#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "CPoint.h"
#include "CVector.h"


CPoint::CPoint(){}
CPoint::CPoint(double x2,double y2,double z2) : x (x2), y (y2), z (z2) {}
CPoint::CPoint(const CPoint &p){
    x = p.x;
    y = p.y;
    z = p.z;
}
double CPoint::getX(){
    return x;
}
double CPoint::getY(){
    return y;
}
double CPoint::getZ(){
    return z;
}
void CPoint::setX(double x2){
    x=x2;
}
void CPoint::setY(double y2){
    y=y2;
}
void CPoint::setZ(double z2){
    z=z2;
}

CPoint CPoint::ProjectOnLine(CPoint pl1,CPoint pl2){
    CVector v = CVector(pl1.getX(), pl1.getY(), pl1.getZ(), pl2.getX(), pl2.getY(), pl2.getZ());
    return ProjectOnLine(v,pl1);
}

CPoint CPoint::ProjectOnLine(CVector v,CPoint pl){

    CVector v3= CVector(x-pl.getX(),y-pl.getY(),z-pl.getZ());//vecteur point a projeter
    
    double n = v.Norme();
    double dist = v.Scalar(v3)/n;
    v.Normalize();

    double xtmp = pl.getX()+v.getX()*dist;
    double ytmp = pl.getY()+v.getY()*dist;
    double ztmp = pl.getZ()+v.getZ()*dist;

    CPoint r = CPoint(xtmp,ytmp,ztmp);
    return r;
}

CPoint CPoint::ProjectOnPlan(CPoint pplan,CVector normalplan){
    CVector v2 = CVector(x, y, z, pplan.getX(), pplan.getY(), pplan.getZ());
    double dist = v2.Scalar(normalplan)/normalplan.Norme();

    double xtmp = x-normalplan.getX()*dist;
    double ytmp = y-normalplan.getY()*dist;
    double ztmp = z-normalplan.getZ()*dist;

    CPoint r = CPoint(xtmp,ytmp,ztmp);
    return r;
}

void CPoint::drawPoint(){
    glBegin(GL_POINTS);
    glVertex3f(x,y,z);
    glEnd();
}
