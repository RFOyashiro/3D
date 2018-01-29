#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "CVector.h"
#include "GL/glut.h" 


class CPoint {
	
	private:
		double x;
		double y;
		double z;

	public:

    CPoint();
    CPoint(double x2,double y2,double z2);
    CPoint(const CPoint &p);

	double getX();
	double getY();
	double getZ();
	void setX(double x2);
	void setY(double y2);
	void setZ(double z2);

    CPoint ProjectOnLine(CPoint pl1,CPoint pl2);
	
    CPoint ProjectOnLine(CVector v,CPoint pl);
    CPoint ProjectOnPlan(CPoint pplan,CVector normalplan);

    void drawPoint();
};
