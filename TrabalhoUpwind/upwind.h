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
    float beta=0;
    float a=0,b=0;
    float threshold=1.0e-14;

    //FUNÇOES
    Upwind(float theta);
    Vector2f ADBQUICKEST(MatrixXf M, int i, int j, Vector2f v,float u_f, float u_g, float theta);
    Vector2f FSLS(MatrixXf M, int i, int j, Vector2f v,float u_f, float u_g);

};

#endif // UPWIND_H
