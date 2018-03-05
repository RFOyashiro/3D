///////////////////////////////////////////////////////////////////////////////
// Imagina
// ----------------------------------------------------------------------------
// IN - Synthèse d'images - Modélisation géométrique
// Auteur : Gilles Gesquière
// ----------------------------------------------------------------------------
// Base du TP 1
// programme permettant de créer des formes de bases.
// La forme représentée ici est un polygone blanc dessiné sur un fond rouge
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>

/* Dans les salles de TP, vous avez généralement accès aux glut dans C:\Dev. Si ce n'est pas le cas, téléchargez les .h .lib ...
Vous pouvez ensuite y faire référence en spécifiant le chemin dans visual. Vous utiliserez alors #include <glut.h>.
Si vous mettez glut dans le répertoire courant, on aura alors #include "glut.h"
*/

#include "GL/freeglut.h"
#include "CVector.h"
#include "CPoint.h"

using namespace std;

// Définition de la taille de la fenêtre
#define WIDTH  480

#define HEIGHT 480

// Définition de la couleur de la fenêtre
#define RED   0
#define GREEN 0
#define BLUE  0
#define ALPHA 1


// Touche echap (Esc) permet de sortir du programme
#define KEY_ESC 27
#define PRECISION 0.1

double nbM = 10;
double nbN = 10;


// Entêtes de fonctions
void init_scene();
void render_scene();
GLvoid initGL();
GLvoid window_display();
GLvoid window_reshape(GLsizei width, GLsizei height);
GLvoid window_key(unsigned char key, int x, int y);


int main(int argc, char **argv)
{
    // initialisation  des paramètres de GLUT en fonction
    // des arguments sur la ligne de commande
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA);

    // définition et création de la fenêtre graphique, ainsi que son titre
    glutInitWindowSize(WIDTH, HEIGHT);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Premier exemple : carré");

    // initialisation de OpenGL et de la scène
    initGL();
    init_scene();

    // choix des procédures de callback pour
    // le tracé graphique
    glutDisplayFunc(&window_display);
    // le redimensionnement de la fenêtre
    glutReshapeFunc(&window_reshape);
    // la gestion des événements clavier
    glutKeyboardFunc(&window_key);

    // la boucle prinicipale de gestion des événements utilisateur
    glutMainLoop();

    return 1;
}

// initialisation du fond de la fenêtre graphique : noir opaque
GLvoid initGL()
{
    glClearColor(RED, GREEN, BLUE, ALPHA);
}

// Initialisation de la scene. Peut servir à stocker des variables de votre programme
// à initialiser
void init_scene()
{
}

// fonction de call-back pour l´affichage dans la fenêtre

GLvoid window_display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    // C'est l'endroit où l'on peut dessiner. On peut aussi faire appel
    // à une fonction (render_scene() ici) qui contient les informations
    // que l'on veut dessiner
    render_scene();

    // trace la scène grapnique qui vient juste d'être définie
    glFlush();
}

// fonction de call-back pour le redimensionnement de la fenêtre

GLvoid window_reshape(GLsizei width, GLsizei height)
{
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // ici, vous verrez pendant le cours sur les projections qu'en modifiant les valeurs, il est
    // possible de changer la taille de l'objet dans la fenêtre. Augmentez ces valeurs si l'objet est
    // de trop grosse taille par rapport à la fenêtre.
    glOrtho(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);


    // toutes les transformations suivantes s´appliquent au modèle de vue
    glMatrixMode(GL_MODELVIEW);
}

// fonction de call-back pour la gestion des événements clavier

GLvoid window_key(unsigned char key, int x, int y)
{
    switch (key) {
        case KEY_ESC:
            exit(1);
            break;

        case '+':
            nbM++;
            nbN++;
            glutPostRedisplay();
            break;

        case '-':
            if (nbM > 4) nbM--;
            if (nbN > 3) nbN--;
            glutPostRedisplay();
            break;

        default:
            printf ("La touche %d n´est pas active.\n", key);
            break;
    }
}

double rand_float(double a, double b) {
    return ((double)rand() / RAND_MAX) * (b - a) + a;
}

CPoint somme(vector<double> asum, vector<CPoint> point){
    double resultx = 0;
    double resulty = 0;
    double resultz = 0;

    CPoint final;

    for (int i = 0;i < asum.size(); i++){
        resultx+=asum[i]*point[i].getX();
        cout << resultx << endl;
        resulty+=asum[i]*point[i].getY();
        resultz+=asum[i]*point[i].getZ();
    }
    final.setX(resultx);
    final.setY(resulty);
    final.setZ(resultz);

    return final;
}

void DrawCurve(vector <CPoint> Points /*,long NbPoint*/) {
    glBegin(GL_LINE_STRIP);
    for (unsigned i (0); i < Points.size(); ++i)
        glVertex3f(Points[i].getX(), Points[i].getY(), Points[i].getZ());
    glEnd();
}

void DrawPoly(vector <CPoint> Points /*,long NbPoint*/) {
    glColor3f(rand_float(0.0, 1.0), rand_float(0.0, 1.0), rand_float(0.0, 1.0));
    glBegin(GL_POLYGON);
    for (unsigned i (0); i < Points.size(); ++i)
        glVertex3f(Points[i].getX(), Points[i].getY(), Points[i].getZ());
    glEnd();
}


vector <CPoint> HermiteCubicCurve(CPoint P0, CPoint P1, CVector V0, CVector V1, long NbU) {
    vector <CPoint> Points;

    for (double i (0); i < NbU; ++i){
        CPoint CurrentPoint;

        double j = (double) i/NbU;

        double F1 = 2 * pow(j, 3) - 3 * pow(j, 2) + 1;
        double F2 = -2 * pow(j, 3) + 3 * pow(j, 2);
        double F3 = pow(j, 3) - 2 * pow(j, 2) + j;
        double F4 = pow(j, 3) - pow(j, 2);

        CurrentPoint.setX(F1 * P0.getX() + F2 * P1.getX() + F3 * V0.getX() + F4 * V1.getX());
        CurrentPoint.setY(F1 * P0.getY() + F2 * P1.getY() + F3 * V0.getY() + F4 * V1.getY());
        CurrentPoint.setZ(F1 * P0.getZ() + F2 * P1.getZ() + F3 * V0.getZ() + F4 * V1.getZ());

        Points.push_back(CurrentPoint);

        cout << CurrentPoint.getX() << " " << CurrentPoint.getY() << " " << CurrentPoint.getZ() << endl;
    }

    return Points;
}

unsigned Factoriel (unsigned value) {
    if (value == 0) return 1;
    return value * Factoriel(value - 1);
}

CPoint SumBern (vector <double> Polynome, vector <CPoint> Point) {
    double ResultX = 0;
    double ResultY = 0;
    double ResultZ = 0;
    for (unsigned i (0); i < Polynome.size(); ++i){
        ResultX += Polynome[i] * Point[i].getX();
        ResultY += Polynome[i] * Point[i].getY();
        ResultZ += Polynome[i] * Point[i].getZ();
    }
    CPoint result;
    result.setX(ResultX);
    result.setY(ResultY);
    result.setZ(ResultZ);

    return result;
}

vector <CPoint> BezierCurveByBernstein(vector <CPoint> Points, long NbU) {
    vector <CPoint> Result;

    for (unsigned i (0); i < NbU; ++i) {
        double j = (double) i / NbU;

        vector <double> PolynomeBernstein;
        for (unsigned k (0); k < Points.size(); ++k) {
            double Facteur1 = (double) Factoriel(Points.size()-1) / (Factoriel(k) * Factoriel(Points.size()-1 - k));
            double Facteur2 = pow(j, k) * pow(1 - j, Points.size()-1 - k);
            PolynomeBernstein.push_back(Facteur1 * Facteur2);
        }
        Result.push_back(SumBern(PolynomeBernstein, Points));
    }

    return Result;
}


vector <vector<CPoint> > SurfaceCylindrique(vector <CPoint> bezier, CVector droite, unsigned nbV, unsigned nbU)
{
    vector <vector<CPoint> > Surface;
    vector<CPoint> Points;

    for(double i (0); i < nbV; i++)
    {
        double v = i / (double) (nbV - 1);

        for(double j (0); j < nbU; j++)
        {
            double u = j / (double) nbU;
            int index = bezier.size() * u;

            double x = bezier[index].getX() + droite.getX() * v;
            double y = bezier[index].getY() + droite.getY() * v;
            double z = bezier[index].getZ() + droite.getZ() * v;

            Points.push_back(CPoint(x,y,z));


            droite.drawLine(bezier[index].getX(), bezier[index].getY(), bezier[index].getZ());
        }

        Surface.push_back(Points);
        Points.clear();

    }
    return Surface;
}

vector <vector<CPoint> > SurfaceReglee(vector <CPoint> Courbe1,  vector <CPoint> Courbe2, unsigned nbV, unsigned nbU)
{
    vector <vector<CPoint> > Surface;
    vector<CPoint> Points;

    for (double i (0); i < nbV; i++) {

        double v = i / (double) nbV ;

        for (double j (0); j < nbU; j++) {
            double u = j / (double) nbU;
            int index = Courbe1.size() * u;

            CVector vect = CVector (Courbe1[i].getX(), Courbe1[i].getY(), Courbe1[i].getZ(), Courbe2[i].getX(), Courbe2[i].getY(), Courbe2[i].getZ());

            double x = (1 - v) * Courbe1[index].getX() + Courbe2[index].getX() * v;
            double y = (1 - v) * Courbe1[index].getY() + Courbe2[index].getY() * v;
            double z = (1 - v) * Courbe1[index].getZ() + Courbe2[index].getZ() * v;

            Points.push_back(CPoint(x, y, z));

            vect.drawLine(Courbe1[i].getX(), Courbe1[i].getY(), Courbe1[i].getZ());
        }

        Surface.push_back(Points);
        Points.clear();

    }
    return Surface;

}

vector<CPoint> BezierCurveByCasteljau(vector<CPoint> points,long nbu)
{
    vector<CPoint> bezierPoints;

    for(double t = 0 ; t < nbu ; t++)
    {
        double k = (double) t/nbu;
        vector<CPoint> tmp1 = points;
        while(tmp1.size()>1)
        {
            vector<CPoint> tmp2;
            for(int i = 0 ; i<tmp1.size()-1 ; i++)
            {
                CPoint Ptmp1 = tmp1[i];
                CPoint Ptmp2 = tmp1[i+1];
                CPoint tmp = CPoint(Ptmp1.getX() * (1-k) + Ptmp2.getX() * k, Ptmp1.getY() * (1-k) + Ptmp2.getY() * k, 0);
                //tmp.drawPoint();
                tmp2.push_back(tmp);
            }
            tmp1 = tmp2;
        }
        bezierPoints.push_back(tmp1[0]);
    }
    return bezierPoints;
}

vector <vector<CPoint> > BezierSurfaceByCasteljau (vector <CPoint> PointsControleU, double nbU, vector <CPoint> PointsControleV, double nbV) {
    vector <vector<CPoint> > Surface;
    vector <CPoint> Points;

    // double pts = (1 - v) * ((1 - u) * P(k, i, j) + u * P(k, i+1, j)) + v * ((1 - u) * P(k, i, j+1) + u * P(k, i+1, j+1));

    if (PointsControleU.size() == PointsControleV.size()) {

        vector <vector <CPoint> > MatriceControle;

        for (unsigned i (0); i < PointsControleV.size(); ++i) {
            vector <CPoint> temp;
            for (unsigned j (0); j < PointsControleU.size(); ++j) {

                if (i == 0) {
                    temp.push_back(PointsControleU[j]);
                }
                else {
                    CPoint P1 = PointsControleV[i - 1];
                    CPoint P2 = PointsControleV[i];



                    CVector vect = CVector(P1.getX(), P1.getY(), P1.getZ(), P2.getX(), P2.getY(), P2.getZ());
                    CPoint ref = MatriceControle[i - 1][j];
                    temp.push_back(CPoint(ref.getX() + vect.getX(), ref.getY() + vect.getY(), ref.getZ() + vect.getZ()));
                }
            }

            MatriceControle.push_back(temp);
        }

        for(double i = 0; i < nbU; i++)
        {
            double v = i / nbU;

            for (double j = 0; j < nbV; ++j)
            {
                double u = j / nbV;

                vector <vector <CPoint> > MatriceTemp = MatriceControle;
                unsigned NbPU = PointsControleU.size();
                unsigned NbPV = PointsControleV.size();

                while (NbPU > 1 && NbPV > 1)
                {


                    vector <CPoint> PointsTemp;
                    vector <vector <CPoint> > MatTemp;

                    for (int k = 0; k < NbPU - 1; k++){
                        for (int l = 0; l < NbPV - 1; l++){

                            CPoint p1 = MatriceTemp[k][l];
                            CPoint p2 = MatriceTemp[k+1][l];
                            CPoint p3 = MatriceTemp[k+1][l+1];
                            CPoint p4 = MatriceTemp[k][l+1];

                            double xa = (1 - u) * p3.getX() + u * p4.getX();
                            double ya = (1 - u) * p3.getY() + u * p4.getY();
                            double za = (1 - u) * p3.getZ() + u * p4.getZ();;

                            double xb = (1 - u) * p2.getX() + u * p1.getX();
                            double yb = (1 - u) * p2.getY() + u * p1.getY();
                            double zb = (1 - u) * p2.getZ() + u * p1.getZ();

                            double xk = (1 - v) * xb + v * xa;
                            double yk = (1 - v) * yb + v * ya;
                            double zk = (1 - v) * zb + v * za;

                            PointsTemp.push_back(CPoint(xk, yk, zk));
                            PointsTemp[PointsTemp.size()-1].drawPoint();
                        }
                        MatTemp.push_back(PointsTemp);
                    }

                    NbPU--;
                    NbPV--;

                    MatriceTemp = MatTemp;

                }
                Points.push_back(MatriceTemp[0][0]);
            }
            Surface.push_back(Points);
            Points.clear();
        }
    }

    return Surface;
}


vector<vector<CPoint> > drawCylindre(){
    double hauteur = 20.0;
    double rayon = 10.0;
    int nbMeridiens = 10;

    vector<CPoint> tmp1;
    vector<CPoint> tmp2;
    vector<vector<CPoint> > result;

    for(double i = 0.0; i < nbMeridiens; i++){
        double omega = 2.0 * M_PI *  i / nbMeridiens;

        double x = rayon * sin(omega);
        double y = rayon * cos(omega);
        double z = 0.0;
        double zb = hauteur;

        CPoint b = CPoint(x, y, z);
        CPoint b2 = CPoint(x, y, zb);

        tmp1.push_back(b);
        tmp2.push_back(b2);


    }
    result.push_back(tmp1);

    result.push_back(tmp2);

    for (int i = 0; i < result[0].size() - 1 ; i++){
        glColor3f(rand_float(0.0, 1.0), rand_float(0.0, 1.0), rand_float(0.0, 1.0));
        glBegin(GL_POLYGON);
        result[0][i].drawPoint();
        result[1][i].drawPoint();
        result[1][i+1].drawPoint();
        result[0][i+1].drawPoint();
        glEnd();

    }

    glBegin(GL_POLYGON);
    result[0][0].drawPoint();
    result[0][result[0].size() - 1].drawPoint();
    result[1][result[0].size() - 1].drawPoint();
    result[1][0].drawPoint();
    glEnd();

    glColor3f(rand_float(0.0, 1.0), rand_float(0.0, 1.0), rand_float(0.0, 1.0));
    glBegin(GL_POLYGON);
    for (int i = 0; i < result[0].size(); i++){
        result[0][i].drawPoint();
    }

    glEnd();
    glBegin(GL_POLYGON);
    result[0][0].drawPoint();
    result[0][result.size()- 1].drawPoint();
    result[0][result.size()- 2].drawPoint();
    glEnd();


    return result;

}

vector<CPoint> drawCone(int n){
    int rayon = 15;
    int hauteur = 20;

    vector<CPoint> result;
    result.push_back(CPoint(0, 0, 20));

    for(double i = 0.0; i < n; i++){
        double omega = 2.0 * M_PI *  i / n;

        double x = rayon * sin(omega);
        double y = rayon * cos(omega);
        double z = 20 + hauteur;

        CPoint b = CPoint(x, y, z);

        result.push_back(b);
    }

    for (int i = 0; i < result.size() - 1 ; i++){
        glColor3f(rand_float(0.0, 1.0), rand_float(0.0, 1.0), rand_float(0.0, 1.0));
        glBegin(GL_POLYGON);
        result[0].drawPoint();
        result[i].drawPoint();
        result[i+1].drawPoint();
        glEnd();
    }

    glColor3f(rand_float(0.0, 1.0), rand_float(0.0, 1.0), rand_float(0.0, 1.0));
    glBegin(GL_POLYGON);
    result[0].drawPoint();
    result[1].drawPoint();
    result[result.size() -1].drawPoint();
    glEnd();


    return result;
}

vector<vector<CPoint> > drawSphere(int Rayon, double nbMeridiens, double nbParallelles){
    vector<vector<CPoint> > Result;

    CPoint Nord = CPoint(0, Rayon, 0);
    CPoint Sud = CPoint(0, -Rayon, 0);
    glBegin(GL_POINTS);
    glColor3f(0, 1, 0);
    glColor3f(1, 0, 0);
    Sud.drawPoint();
    glColor3f(0, 0, 1);
    glEnd();

    for (int i = 0; i < nbParallelles; i++){
        double phi =  M_PI * i / nbParallelles;
        vector<CPoint> TempPoints;

        for (int j = 0; j < nbMeridiens; j++){
            double theta = 2 * M_PI * j / nbMeridiens;
            double x = Rayon * sin(phi) * cos(theta);
            double z = Rayon * sin(phi) * sin(theta);
            double y = Rayon * cos(phi);
            CPoint t = CPoint(x, y, z);
            glBegin(GL_POINTS);
            t.drawPoint();
            glEnd();
            TempPoints.push_back(t);
        }

        Result.push_back(TempPoints);
        TempPoints.clear();
    }

    for (int i = 0; i < Result.size() - 1; i++){

        for (int j = 0; j < Result[i].size() - 1; j++){
            glColor3f(rand_float(0.0, 1.0), rand_float(0.0, 0.5), rand_float(0.0, 0.5));
            glBegin(GL_POLYGON);
            Result[i][j].drawPoint();
            Result[i+1][j].drawPoint();
            Result[i+1][j+1].drawPoint();
            Result[i][j+1].drawPoint();
            glEnd();
        }

        glColor3f(rand_float(0.0, 0.5), rand_float(0.0, 1.0), rand_float(0.0, 0.5));
        glBegin(GL_POLYGON);
        Result[i][0].drawPoint();
        Result[i+1][0].drawPoint();
        Result[i+1][nbM-1].drawPoint();
        Result[i][nbM-1].drawPoint();
        glEnd();
    }

    for (int i = 0; i < nbMeridiens; i++){
        glColor3f(rand_float(0.0, 0.5), rand_float(0.0, 0.5), rand_float(0.0, 1.0));
        glBegin(GL_POLYGON);
        Nord.drawPoint();
        Result[0][i].drawPoint();
        Result[0][i+1].drawPoint();
        glEnd();
    }

    glBegin(GL_POLYGON);
    Nord.drawPoint();
    Result[0][0].drawPoint();
    Result[0][nbParallelles-1].drawPoint();
    glEnd();

    for (int i = 0; i < nbMeridiens; i++){
        glColor3f(rand_float(0.0, 1.0), rand_float(0.0, 1.0), rand_float(0.0, 1.0));
        glBegin(GL_POLYGON);
        Sud.drawPoint();
        Result[nbParallelles-1][i].drawPoint();
        Result[nbParallelles-1][i+1].drawPoint();
        glEnd();
    }

    glBegin(GL_POLYGON);
    Sud.drawPoint();
    Result[nbParallelles-1][0].drawPoint();
    Result[nbParallelles-1][nbMeridiens-1].drawPoint();
    glEnd();

    return Result;
}

void DisplayVoxel (CPoint Centre, double length)
{
    vector <CPoint> VoxelP1;

    VoxelP1.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));
    VoxelP1.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() - length / 2, Centre.getZ() + length / 2));
    VoxelP1.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() - length / 2, Centre.getZ() + length / 2));
    VoxelP1.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));
    VoxelP1.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));

    DrawPoly(VoxelP1);

    vector <CPoint> VoxelP2;

    VoxelP2.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));
    VoxelP2.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() - length / 2, Centre.getZ() - length / 2));
    VoxelP2.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() - length / 2, Centre.getZ() - length / 2));
    VoxelP2.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));
    VoxelP2.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));

    DrawPoly(VoxelP2);

    for (unsigned i (0); i < VoxelP1.size() - 1; ++i) {
        glBegin(GL_POLYGON);
        glVertex3f(VoxelP1[i].getX(), VoxelP1[i].getY(), VoxelP1[i].getZ());
        glVertex3f(VoxelP2[i].getX(), VoxelP2[i].getY(), VoxelP2[i].getZ());
        glVertex3f(VoxelP2[i+1].getX(), VoxelP2[i+1].getY(), VoxelP2[i+1].getZ());
        glVertex3f(VoxelP1[i+1].getX(), VoxelP1[i+1].getY(), VoxelP1[i+1].getZ());
        glEnd();
    }
}

void getVoxels (CPoint Point, double length) {
    vector <CPoint> Voxel;

    Voxel.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));
    Voxel.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() - length / 2, Centre.getZ() + length / 2));
    Voxel.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() - length / 2, Centre.getZ() + length / 2));
    Voxel.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));
    Voxel.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() + length / 2));
    Voxel.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));
    Voxel.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() - length / 2, Centre.getZ() - length / 2));
    Voxel.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() - length / 2, Centre.getZ() - length / 2));
    Voxel.push_back(CPoint(Centre.getX() - length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));
    Voxel.push_back(CPoint(Centre.getX() + length / 2, Centre.getY() + length / 2, Centre.getZ() - length / 2));

    return Voxel;
}

vector<CPoint> BreakVoxel(vecto<vector<CPoint>> & Voxel, double Resolution, unsigned index) {
    //calculer nouveau centre => voir feuille TP
    //taille = resolution / 2


}

void DisplaySphereVolumic (CPoint Centre, double Rayon, double Resolution)
{
    vector<vector<CPoint> > Result;

    vector<vector<CPoint> > Voxels;

    for (int i = 0; i < nbN; i++){
        double phi =  M_PI * i / nbN;
        vector<CPoint> TempPoints;

        for (int j = 0; j < nbM; j++){
            double theta = 2 * M_PI * j / nbM;
            double x = Rayon * sin(phi) * cos(theta);
            double z = Rayon * sin(phi) * sin(theta);
            double y = Rayon * cos(phi);
            CPoint t = CPoint(x + Centre.getX(), y + Centre.getY(), z + + Centre.getZ());
            TempPoints.push_back(t);
        }

        Result.push_back(TempPoints);
        TempPoints.clear();
    }

    for (unsigned i (0); i < Result.size(); ++i)
        for (unsigned j (0); j < Result[i].size(); ++j)
            Voxels.push_back(getVoxels(Result[i][j], Resolution));


    for (unsigned i (0); i < Voxels.size(); ++i) {
        bool isOut = false;
        bool isIn = true;
        for (unsigned j (0); j < Voxels[i].size(); ++j) {
            CVector tmp = CVector (Voxels[i][j].getX(), Voxels[i][j].getY(), Voxels[i][j].getZ(), Centre.getX(), Centre.getY(), Centre.getZ());

            if (tmp.Norme() > Rayon) {
                isOut = true;
                isIn = false;
                continue;
            }

            if (isOut && tmp.Norme() < Rayon) {
                BreakVoxel(Voxels, Resolution, i);
                Voxels.erase(i);
                isOut = false;
                break;
            }

        }
        if (isOut) Voxels.erase(i);
        if (isIn) DrawPoly(Voxels[i]);
    }


}

void DisplayCylinderVolumic (CPoint Origin, CVector Vect, double Rayon, double Resolution)
{
    vector<vector<CPoint> > Result;

    for (unsigned i (0); i < Vect.Norme(); ++i) {
            vector <CPoint> tmp1;

            for(double j = 0.0; j < nbM; j++){
                double omega = 2.0 * M_PI *  j / nbM;

                double x = Rayon * sin(omega);
                double y = Rayon * cos(omega);
                double z = i;

                CPoint b = CPoint(x + Origin.getX(), y + Origin.getY(), z + Origin.getZ());

                tmp1.push_back(b);
            }
            Result.push_back(tmp1);
            tmp1.clear();
    }

    for (unsigned i (0); i < Result.size(); ++i)
        for (unsigned j (0); j < Result[i].size(); ++j)
            DisplayVoxel(Result[i][j], Resolution);

}

void Display_INTERSECT_CS (CPoint CentreS, double RayonS, CPoint OriginC, CVector AxeC, double RayonC, double Resolution) {
    vector<vector<CPoint> > Sphere;
    vector<vector<CPoint> > Cylindre;

/*    Sphere = PointSphere(CentreS, RayonS);
    Cylindre = PointCylindre(OriginC, AxeC, RayonC);
*/





}



//////////////////////////////////////////////////////////////////////////////////////////
// Fonction que vous allez modifier afin de dessiner
/////////////////////////////////////////////////////////////////////////////////////////
void render_scene()
{
    glOrtho(-51.0, 51.0, -51.0, 51.0, -51.0, 51.0);
    gluLookAt(15.0, 20.0, 15.0, 0, 10, 0, 0, 1, 0);

    nbN = 20;
    nbM = 20;

    CPoint Centre = CPoint(50, 0, 0);
    CPoint Centre2 = CPoint(0, 0, 0);
    CVector Vect = CVector(0, 0, 0, 20, 20, 20);

    DisplaySphereVolumic(Centre, 30, 3);
    DisplaySphereVolumic(Centre2, 30, 3);

}























































