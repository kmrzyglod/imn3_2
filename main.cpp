#include <iostream>
#include "imnmath.h"
#include "Point.h"
#include "FlagMatrix.h"
#include "Relaxation.h"
#include "RelaxationNeumann.h"

using namespace std;

#define XSIZE 200
#define YSIZE 100
#define TOL 1e-6

double boundaryConditionDirichlet(int x, int y ) {
    x++;
    y++;
    double u0 = 1;
    int  xmin, ymin;
    xmin = ymin = 1;
    int xsize  = 200;
    int ysize = 100;
    if((x == xmin || x == xsize) && y>ymin && y<ysize) {
        return u0*y;
    }
    else if (y <= 30) {
        return u0*ymin;
    }
    else {
        return  u0*ysize;
    }
}
double boundaryConditionNeumann(int x, int y ) {
    x++;
    y++;
    double u0 = 1;
    int  xmin, ymin;
    xmin = ymin = 1;
    int xsize  = 200;
    int ysize = 100;
    if((x == xmin || x == xsize) && y>ymin && y<ysize) {
        return u0*x;
    }
    else if (y <= 30) {
        return u0*ymin;
    }
    else {
        return  u0*ysize;
    }
}

void zad1() {
    vector<Point> obstacle1 = {{85,90}, {100,90}, {100,70}, {115,70}, {115, 100}, {85, 100}, {85,90}};
    vector<Point> obstacle2 = {{85,1}, {85,10}, {100,10}, {100,30}, {115, 30}, {115, 1}, {85,1}};
    FlagMatrix flagMatrix = FlagMatrix(XSIZE, YSIZE);
    flagMatrix.DrawObstacle(obstacle1);
    flagMatrix.DrawObstacle(obstacle2);
    flagMatrix.SetBorders();
    Relaxation relaxation = Relaxation(XSIZE,YSIZE,&flagMatrix, boundaryConditionDirichlet);
    relaxation.NextIteration();
    relaxation.NextIteration();
    while(relaxation.GetTolerance()>=TOL) {
        relaxation.NextIteration();
    }
    cout  << relaxation.GetIteration();
    relaxation.SaveResults();

}
void zad2() {
    vector<Point> obstacle1 = {{85,90}, {100,90}, {100,70}, {115,70}, {115, 100}, {85, 100}, {85,90}};
    vector<Point> obstacle2 = {{85,1}, {85,10}, {100,10}, {100,30}, {115, 30}, {115, 1}, {85,1}};
    FlagMatrix flagMatrix = FlagMatrix(XSIZE, YSIZE);
    flagMatrix.DrawObstacle(obstacle1);
    flagMatrix.DrawObstacle(obstacle2);
    flagMatrix.SetBorders();
    flagMatrix.ConvertToNeumann();
    RelaxationNeumann relaxation = RelaxationNeumann(XSIZE,YSIZE,&flagMatrix, boundaryConditionNeumann);
    relaxation.NextIteration();
    relaxation.NextIteration();
    while(relaxation.GetTolerance()>=TOL) {
        relaxation.NextIteration();
    }
    cout  << relaxation.GetIteration();
    relaxation.SaveResults();

}



int main() {
    zad1();
    zad2();
    return 0;
}