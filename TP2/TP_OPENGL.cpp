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

#include "GL/glut.h" 
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
    glutCreateWindow("Courbe Cubique d'Hermite");

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
    glOrtho(-4.0, 4.0, -4.0, 4.0, -4.0, 4.0);

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

        default:
            printf ("La touche %d n'est pas active.\n", key);
            break;
    }
}

void DrawCurve(vector <CPoint> Points /*,long NbPoint*/) {
    glBegin(GL_LINE_STRIP);
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

vector <CPoint> BezierCurveByCasteljau(vector <CPoint> Points, long NbU) {
    vector <CPoint> Result;

    CPoint Point;

    for (unsigned i (0); i < NbU; ++i) {
        double j = (double) i / NbU;

        vector <CPoint> TempPoint = Points;

        while (TempPoint.size() > 1) {

            vector <CPoint> NewPoint;

            for (unsigned k (0); k < TempPoint.size() - 1; ++k) {
                CPoint Point;
                Point.setX((1 - j) * TempPoint[k].getX() + j * TempPoint[k + 1].getX());
                Point.setY((1 - j) * TempPoint[k].getY() + j * TempPoint[k + 1].getY());
                Point.setZ((1 - j) * TempPoint[k].getZ() + j * TempPoint[k + 1].getZ());
                NewPoint.push_back(Point);
            }
            TempPoint.clear();
            TempPoint = NewPoint;
        }
        Point.setX(j * TempPoint[0].getX());
        Point.setY(j * TempPoint[0].getY());
        Point.setZ(j * TempPoint[0].getZ());
        Result.push_back(Point);
    }
    return Result;

}


//////////////////////////////////////////////////////////////////////////////////////////
// Fonction que vous allez modifier afin de dessiner
/////////////////////////////////////////////////////////////////////////////////////////
void render_scene()
{
    unsigned NbU = 30;

    CPoint P0 = CPoint(-1, -1, 0);
    CPoint P1 = CPoint(0, 2, 0);
    CPoint P2 = CPoint(1, 3, 0);
    CPoint P3 = CPoint(2, 2, 0);
    CPoint P4 = CPoint(3, -1, 0);

    vector <CPoint> PointsControle;
    PointsControle.push_back(P0);
    PointsControle.push_back(P1);
    PointsControle.push_back(P2);
    PointsControle.push_back(P3);
    PointsControle.push_back(P4);

    vector <CPoint> Points = BezierCurveByCasteljau(PointsControle, NbU);

    glColor3f(1.0, 0.0, 0.0); //red
    DrawCurve(PointsControle);

    glColor3f(1.0, 1.0, 1.0); //white
    DrawCurve(Points);
}
