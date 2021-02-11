#ifndef UPWIND_H
#define UPWIND_H
#include "../Eigen/Eigen/dense"
#include <iostream>

using namespace Eigen;
using namespace std;

class Upwind
{
private:
public:
    //VARIAVEIS
    int tipoDeUpwind=0;
    double beta=0;
    double a=0,b=0;
    double threshold=1.0e-14;

    //FUNÇOES
    Upwind(double theta);
    Vector2f ADBQUICKEST(MatrixXf M, int i, int j, Vector2f v,double u_f, double u_g, double theta);
    Vector2f FSLS(MatrixXf M, int i, int j, Vector2f v,double u_f, double u_g);

};

#endif // UPWIND_H
