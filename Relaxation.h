//
// Created by Kamil on 08.05.2016.
//
#ifndef IMN3_2_RELAXATION_H
#define IMN3_2_RELAXATION_H

#include <c++/functional>
#include "FlagMatrix.h"

class Relaxation {
private:
    double** _relaxationMatrix;
    int _xsize, _ysize, _iter;
    double _nowIntegral, _prevIntegral, _tolerance;
    FlagMatrix* _flagMatrix;
    function < double(int x, int y) > _boundaryConditionFn;


    void makeRelaxation(){
        for(int i=0;i<_xsize;i++) {
            for(int j=0;j<_ysize;j++) {
                if(_flagMatrix->GetMatrix()[i][j] == 0)
                {
                    _relaxationMatrix[i][j] = (_relaxationMatrix[i-1][j] + _relaxationMatrix[i][j-1]
                                              + _relaxationMatrix[i+1][j] + _relaxationMatrix[i][j+1])/4.0;
                }
            }
        }
    }

    void calculateIntegral() {
        double integral = 0;
        double a,b;
        for(int i=0;i<_xsize;i++) {
            for(int j=0;j<_ysize;j++) {
                if(_flagMatrix->GetMatrix()[i][j] == 0)
                {
                    a = pow((_relaxationMatrix[i+1][j] - _relaxationMatrix[i-1][j]), 2);
                    b = pow((_relaxationMatrix[i][j+1] - _relaxationMatrix[i][j-1]), 2);
                    integral+= a + b;
                }
            }
        }
        _nowIntegral = 1.0/8.0 * integral;
    }
    void calculateTolerance() {
        _tolerance = fabs((_nowIntegral - _prevIntegral)/_prevIntegral);
    }
    void setBoundaryConditions() {
        for(int i=0;i<_xsize;i++) {
            for(int j=0;j<_ysize;j++) {
                if(_flagMatrix->GetMatrix()[i][j] != 0)
                {
                    _relaxationMatrix[i][j] = _boundaryConditionFn(i, j);
                }
            }
        }
    }
public:
    Relaxation(int xsize, int ysize, FlagMatrix* flagMatrix,  function < double( int x, int y) > boundaryConditionFn):
            _xsize(xsize), _ysize(ysize), _flagMatrix(flagMatrix), _boundaryConditionFn(boundaryConditionFn) {
        _relaxationMatrix = imnd::matrix(_xsize, _ysize);
        Reset();
        //PrintMatrixToFile("test.txt");
        //SaveFlagsMatrixToPNGFile("test_conditions.png");
    }
    void Reset() {
        imnd::set_matrix(_relaxationMatrix, _xsize, _ysize, 0);
        _nowIntegral = 0;
        _prevIntegral = 0;
        _iter = 0;
        setBoundaryConditions();
    }
    void NextIteration() {
        makeRelaxation();
        calculateIntegral();
        if(_iter > 0) {
            calculateTolerance();
        }
        _prevIntegral = _nowIntegral;
        _iter++;
    }
    int GetIteration() {
        return _iter;
    }
    double GetTolerance() {
        return _tolerance;
    }
    void SaveFlagsMatrixToPNGFile(const char *fileName) {
        imnd::plot_2d_system(fileName, _relaxationMatrix, _xsize, _ysize,  1, 1);
    }
    void SaveResults() {
        imnd::plot_params.stype = GNUPLOT_CONTOUR | GNUPLOT_PM3D;
        SaveFlagsMatrixToPNGFile("test_relaxation.png");
    }
    void PrintMatrixToFile(const char* fileName) {
        ofstream oFile;
        oFile.open (fileName);
        for(int i=0;i<_xsize;i++) {
            for(int j=0;j<_ysize;j++) {
                oFile << _relaxationMatrix[i][j] << ", ";
            }
            oFile << "\b\b\n";
        }
        oFile.close();
    }

};

#endif //IMN3_2_RELAXATION_H
