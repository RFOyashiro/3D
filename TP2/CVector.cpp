#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "CVector.h"
#include "GL/glut.h"

CVector::CVector(){}
CVector::CVector(double x2,double y2,double z2) : x (x2), y (y2), z (z2) {}
CVector::CVector(const CVector & p) {
   x = p.x;
   y = p.y;
   z = p.z;
}

CVector::CVector(double x1, double y1, double z1, double x2, double y2, double z2){
    x=x1-x2;
    y=y1-y2;
    z=z1-z2;
}

double CVector::getX(){
    return x;
}
double CVector::getY(){
    return y;
}
double CVector::getZ(){
    return z;
}
void CVector::setX(double x2){
        x=x2;
}
void CVector::setY(double y2){
        y=y2;
}
void CVector::setZ(double z2){
        z=z2;
}
double CVector::Norme(){
    return sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0));
}
void CVector::Normalize(){
    double n=Norme();
    x=x/n;
    y=y/n;
    z=z/n;
}
double CVector::Scalar(CVector v2){
    return x*v2.getX()+y*v2.getY()+z*v2.getZ();
}

CVector CVector::Vectoriel(CVector v2){
    double xtmp = y*v2.getZ()-z*v2.getY();
    double ytmp = z*v2.getX()-x*v2.getZ();
    double ztmp = x*v2.getY()-y*v2.getX();

        CVector tmp = CVector(xtmp,ytmp,ztmp);

	return tmp;
}

double CVector::Angle(CVector v2){
    return acos(this->Scalar(v2)/(Norme()*v2.Norme()));
}

void CVector::drawLine(double x1, double y1, double z1){
	glBegin(GL_LINES);
        glVertex3f(x1,y1,z1);
        glVertex3f(x1+x,y1+y,z1+z);
	glEnd();
}
