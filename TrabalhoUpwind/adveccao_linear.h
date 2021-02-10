#ifndef ADVECCAO_LINEAR_H
#define ADVECCAO_LINEAR_H
#include "../Eigen/Eigen/dense"
#include <iostream>

using namespace Eigen;
using namespace std;


class Adveccao_Linear
{
public:
    float visc,
          deltaX,
          delta_t,
          comprimento,
          tempo,
          contornoInicio,
          contornoFim,
          theta=0;
    int contX;
    int cont_t;

    Adveccao_Linear();
    MatrixXf CalculaEquacao(VectorXf X);
    MatrixXf InsereCondicoesIniciais(MatrixXf M, VectorXf X);
};

#endif // ADVECCAO_LINEAR_H
