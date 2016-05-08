#include <iostream>
#include "imnmath.h"
#include "Point.h"
#include "FlagMatrix.h"
#include "Relaxation.h"

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
    else if (y == ymin) {
        return u0*ymin;
    }
    else {
        return  u0*ysize;
    }
}

void zad1() {
    vector<Point> obstacle = {{85,85}, {100,85}, {100,70}, {115,70}, {115, 100}, {85, 100}, {85,85}};
    FlagMatrix flagMatrix = FlagMatrix(XSIZE, YSIZE);
    flagMatrix.DrawObstacle(obstacle);
    Relaxation relaxation = Relaxation(XSIZE,YSIZE,&flagMatrix, boundaryConditionDirichlet);
    relaxation.NextIteration();
    relaxation.NextIteration();
    //cout << relaxation.GetTolerance() << "\n";
    while(relaxation.GetTolerance()>=TOL) {
       // cout << relaxation.GetTolerance() << "\n";
        relaxation.NextIteration();
    }
    cout  << relaxation.GetIteration();
    relaxation.SaveResults();

}
void zad2() {

}



int main() {
    zad1();
    return 0;
}