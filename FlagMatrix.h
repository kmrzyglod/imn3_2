//
// Created by Kamil on 08.05.2016.
//
#ifndef IMN3_FLAGMATRIX_H
#define IMN3_FLAGMATRIX_H

#include <vector>
#include <math.h>
#include <fstream>
#include "imnmath.h"

typedef imn<double> imnd;
typedef imn<int> imni;

using namespace std;
class FlagMatrix {
private:
    int** _flagMatrix;
    int _xsize, _ysize;

    void fillObstacle() {
        bool insideSentry = false;
        for(int i=0;i<_xsize;i++) {
            insideSentry = false;
            for(int j=0;j<_ysize-1;j++) {
                if(_flagMatrix[i][j] == 1) {
                    if(insideSentry) {
                        insideSentry = false;
                    }
                    else if(_flagMatrix[i][j+1] == 0 && (_flagMatrix[i][j-1] == 0 || j==0 )) {
                        insideSentry = true;
                    }
                    else {
                        continue;
                    }
                } else {
                    if(insideSentry) {
                        _flagMatrix[i][j] = 1;
                    }
                }
            }
        }
    }
    //Bresenham's Line algorithm implementation
    void bhmLine(int x1, int y1, int x2, int y2, int c)
    {
        int x,y,dx,dy,dx1,dy1,px,py,xe,ye,i;
        dx=x2-x1;
        dy=y2-y1;
        dx1=fabs(dx);
        dy1=fabs(dy);
        px=2*dy1-dx1;
        py=2*dx1-dy1;
        if(dy1<=dx1)
        {
            if(dx>=0)
            {
                x=x1;
                y=y1;
                xe=x2;
            }
            else
            {
                x=x2;
                y=y2;
                xe=x1;
            }
            _flagMatrix[x-1][y-1] = 1;
            for(i=0;x<xe;i++)
            {
                x=x+1;
                if(px<0)
                {
                    px=px+2*dy1;
                }
                else
                {
                    if((dx<0 && dy<0) || (dx>0 && dy>0))
                    {
                        y=y+1;
                    }
                    else
                    {
                        y=y-1;
                    }
                    px=px+2*(dy1-dx1);
                }
                _flagMatrix[x-1][y-1] = 1;
            }
        }
        else
        {
            if(dy>=0)
            {
                x=x1;
                y=y1;
                ye=y2;
            }
            else
            {
                x=x2;
                y=y2;
                ye=y1;
            }
            _flagMatrix[x-1][y-1] = 1;
            for(i=0;y<ye;i++)
            {
                y=y+1;
                if(py<=0)
                {
                    py=py+2*dx1;
                }
                else
                {
                    if((dx<0 && dy<0) || (dx>0 && dy>0))
                    {
                        x=x+1;
                    }
                    else
                    {
                        x=x-1;
                    }
                    py=py+2*(dx1-dy1);
                }
                _flagMatrix[x-1][y-1] = 1;
            }
        }
    }
public:
    FlagMatrix(int xsize, int ysize): _xsize(xsize), _ysize(ysize){
        _flagMatrix = imn<int>::matrix(_xsize, _ysize);
        imni::set_matrix(_flagMatrix, _xsize, _ysize, 0);
    }

    void SetBorders() {
        imni::set_matrix_col(_flagMatrix, 0, _xsize, 1);
        imni::set_matrix_col(_flagMatrix, _ysize-1, _xsize, 1);
        imni::set_matrix_row(_flagMatrix, 0, _ysize, 1);
        imni::set_matrix_row(_flagMatrix, _xsize-1, _ysize, 1);

    }

    void ConvertToNeumann() {
        for(int i=1;i<_xsize-1;i++) {
            for(int j=0;j<_ysize;j++) {
                _flagMatrix[i][j]*=2;
            }
        }
        _flagMatrix[0][_ysize-1]*=2;
        _flagMatrix[_xsize-1][_ysize-1]*=2;
        _flagMatrix[_xsize-1][0]*=2;
        _flagMatrix[0][0]*=2;
    }

    void DrawObstacle(vector<Point>&vertices){
        for(int i=0;i<vertices.size()-1;i++){
            bhmLine(vertices[i].GetX(), vertices[i].GetY(), vertices[i+1].GetX(), vertices[i+1].GetY(), 1);
        }
        fillObstacle();
    }

    void SaveFlagsMatrixToPNGFile(const char *fileName) {
        imni::plot_2d_system(fileName, _flagMatrix, _xsize, _ysize,  1, 1);
    }

    void PrintMatrixToFile(const char* fileName) {
        ofstream oFile;
        oFile.open (fileName);
        for(int i=0;i<_xsize;i++) {
            for(int j=0;j<_ysize;j++) {
                oFile << _flagMatrix[i][j] << ", ";
            }
            oFile << "\b\b\n";
        }
        oFile.close();
    }

    int** GetMatrix() {
        return _flagMatrix;
    }

    ~FlagMatrix() {
        imni::free_matrix(_flagMatrix, _xsize);
    }

};

#endif //IMN3_FLAGMATRIX_H
