//
// Created by Kamil on 08.05.2016.
//

#ifndef IMN3_2_RELAXATIONNEUMANN_H
#define IMN3_2_RELAXATIONNEUMANN_H
#include <c++/functional>
#include "FlagMatrix.h"

class RelaxationNeumann: public Relaxation {
protected:
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

public:
    RelaxationNeumann(int xsize, int ysize, FlagMatrix* flagMatrix,  function < double( int x, int y) > boundaryConditionFn):  Relaxation(xsize, ysize, flagMatrix, boundaryConditionFn) {}

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

};
#endif //IMN3_2_RELAXATIONNEUMANN_H
