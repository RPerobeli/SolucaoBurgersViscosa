#ifndef BURGERS_VISCOSA_H
#define BURGERS_VISCOSA_H
#include "../Eigen/Eigen/dense"
#include <iostream>

using namespace Eigen;
using namespace std;


class Burgers_Viscosa
{
  private:
  public:
    //VARIAVEIS
    double visc,
          deltaX,
          delta_t,
          comprimento,
          tempo,
          contornoInicio,
          contornoFim,
          theta=0;


    int contX;
    int cont_t;

    Burgers_Viscosa();
    MatrixXf CalculaEquacao(VectorXf X);
    MatrixXf InsereCondicoesIniciais(MatrixXf M, VectorXf X);



};

#endif // BURGERS_VISCOSA_H
