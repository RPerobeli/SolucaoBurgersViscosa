#include "adveccao_linear.h"
#include "upwind.h"

#define pi 3.141589

using namespace Eigen;
using namespace std;

Adveccao_Linear::Adveccao_Linear()
{
    cout<<"viscosidade: "<<endl;
    cin >> this->visc;
    cout<<"Theta: "<<endl;
    cin >> this->theta;
    //cout<<"Comprimento X total: "<<endl;
    //cin >> this->comprimento;
    this->comprimento = 1;
    //cout<<"Passo em X: "<<endl;
    //cin >> this->deltaX;
    this->deltaX=0.01;
    cout<<"Tempo total: "<<endl;
    cin >> this->tempo;
    //this->tempo = 0.5
    //cout<<"Passo do tempo: "<<endl;
    //cin >> this->delta_t;
    this->delta_t=theta*deltaX;

}
