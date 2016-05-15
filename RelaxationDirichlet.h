//
// Created by Kamil on 15.05.2016.
//

#include "Relaxation.h"

#ifndef IMN3_2_RELAXATIONDIRICHLET_H
#define IMN3_2_RELAXATIONDIRICHLET_H

class RelaxationDirichlet: public Relaxation {
protected:
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

public:
    RelaxationDirichlet(int xsize, int ysize, FlagMatrix* flagMatrix,  function < double( int x, int y) > boundaryConditionFn): Relaxation(xsize, ysize, flagMatrix, boundaryConditionFn) {}

    void SaveResults() {
        imnd::plot_params.title = "Zad1. Przeplyw potencjalny";
        imnd::plot_params.stype = GNUPLOT_CONTOUR | GNUPLOT_PM3D;
        SaveFlagsMatrixToPNGFile("RelaxationZad1.png");
        ofstream oFile;
        oFile.open ("IntegralZad1.txt");
        for(IntegralGraphPoint el:  _integralGraph) {
            oFile << el.GetX() << " " << el.GetY() << "\n";
        }
        oFile.close();
        system("gnuplot integralZad1.plt");
    }

};

#endif //IMN3_2_RELAXATIONDIRICHLET_H
