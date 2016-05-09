//
// Created by Kamil on 08.05.2016.
//

#ifndef IMN3_2_RELAXATIONNEUMANN_H
#define IMN3_2_RELAXATIONNEUMANN_H
#include <c++/functional>
#include "FlagMatrix.h"

class RelaxationNeumann {
private:
    double** _relaxationMatrix;
    int _xsize, _ysize, _iter;
    double _nowIntegral, _prevIntegral, _tolerance;
    FlagMatrix* _flagMatrix;
    function < double(int x, int y) > _boundaryConditionFn;
    vector<IntegralGraphPoint> _integralGraph;

    void makeRelaxation(){
        //vector<double> neighboors;
        double neighbourPotential = 0;
        int neighbourConuter = 0;
        for(int i=0;i<_xsize;i++) {
            for(int j=0;j<_ysize;j++) {
                if(_flagMatrix->GetMatrix()[i][j] == 0)
                {

                        _relaxationMatrix[i][j] = (_relaxationMatrix[i - 1][j] + _relaxationMatrix[i][j - 1]
                                                   + _relaxationMatrix[i + 1][j] + _relaxationMatrix[i][j + 1]) / 4.0;
                }
                else if(_flagMatrix->GetMatrix()[i][j] == 2){
                    neighbourConuter = 0;
                    neighbourPotential = 0;
                    if((i+1) < _xsize &&  _flagMatrix->GetMatrix()[i+1][j] == 0 )
                    {
                        neighbourConuter++;
                        neighbourPotential+=_relaxationMatrix[i+1][j];
                    }
                    if((i-1) >= 0 &&  _flagMatrix->GetMatrix()[i-1][j] == 0 )
                    {
                        neighbourConuter++;
                        neighbourPotential+=_relaxationMatrix[i-1][j];
                    }
                    if((j+1) < _ysize &&  _flagMatrix->GetMatrix()[i][j+1] == 0 )
                    {
                        neighbourConuter++;
                        neighbourPotential+=_relaxationMatrix[i][j+1];
                    }
                    if((j-1) >= 0 &&  _flagMatrix->GetMatrix()[i][j-1] == 0 )
                    {
                        neighbourConuter++;
                        neighbourPotential+=_relaxationMatrix[i][j-1];
                    }
                    if(neighbourConuter > 0) {
                        _relaxationMatrix[i][j] = neighbourPotential/neighbourConuter;
                    } else {
                       // _relaxationMatrix[i][j] = 0;
                    }
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
                if(_flagMatrix->GetMatrix()[i][j] == 1)
                {
                    _relaxationMatrix[i][j] = _boundaryConditionFn(i, j);
                }
            }
        }
    }
public:
    RelaxationNeumann(int xsize, int ysize, FlagMatrix* flagMatrix,  function < double( int x, int y) > boundaryConditionFn):
            _xsize(xsize), _ysize(ysize), _flagMatrix(flagMatrix), _boundaryConditionFn(boundaryConditionFn) {
        _relaxationMatrix = imnd::matrix(_xsize, _ysize);
        Reset();
        //PrintMatrixToFile("test3.txt");
        //SaveFlagsMatrixToPNGFile("test3.png");

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
            _integralGraph.push_back({_iter, _nowIntegral});

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
        imnd::plot_params.title = "Zad2. Przeplyw potencjalny";
        imnd::plot_params.stype = GNUPLOT_CONTOUR | GNUPLOT_PM3D;
        SaveFlagsMatrixToPNGFile("RelaxationZad2.png");
        ofstream oFile;
        oFile.open ("IntegralZad2.txt");
        for(IntegralGraphPoint el:  _integralGraph) {
            oFile << el.GetX() << " " << el.GetY() << "\n";
        }
        oFile.close();
        system("gnuplot integralZad2.plt");
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
#endif //IMN3_2_RELAXATIONNEUMANN_H
