#include "burgers_viscosa.h"
#include "upwind.h"

#define pi 3.141589

using namespace Eigen;
using namespace std;

Burgers_Viscosa::Burgers_Viscosa()
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

//Burgers_Viscosa::~Burgers_Viscosa(){}

MatrixXf Burgers_Viscosa::CalculaEquacao(VectorXf X)
{

    MatrixXf U(cont_t,contX); //matriz com as Velocidades

    U = this->InsereCondicoesIniciais(U, X);

    //começando o calculo

    float uBarra_f = 0, uBarra_g = 0;

    Upwind upwind(theta);
    Vector2f u_gf = {0,0};

    switch(upwind.tipoDeUpwind)
    {
        case(0):
        {
            cout<<"erro"<<endl;
            break;
        }
        case(1):
        {
            cout<<"FSLS"<<endl;
            for(int n=0; n<U.rows()-1;n++)
            {
                for(int i=1; i<U.cols()-1;i++)
                {
                    //loop do espaço
                    uBarra_f = (U(n,i+1)+U(n,i))*0.5;
                    uBarra_g = (U(n,i)+U(n,i-1))*0.5;

                    u_gf = upwind.FSLS(U,n,i, u_gf,uBarra_f,uBarra_g);
                    float p1 = 0.5*delta_t*(1.0/deltaX)*(uBarra_f* u_gf(1) - uBarra_g * u_gf(0));
                    float p2 = (delta_t*visc/ (deltaX*deltaX))*(U(n,i+1)-2.0*U(n,i) + U(n,i-1));
                    U(n+1,i)= U(n,i) - p1 +  p2;

                }
            }
            cout << "deu certo fsls" <<endl;
            break;
        }
        case(2):
        {
            cout << "ADBQUICKEST"<<endl;
            for(int n=0; n<U.rows()-1;n++)
            {
                for(int i=1; i<U.cols()-1;i++)
                {
                    //loop do espaço
                    uBarra_f = (U(n,i+1)+U(n,i))*0.5;
                    uBarra_g = (U(n,i)+U(n,i-1))*0.5;

                    u_gf = upwind.ADBQUICKEST(U,n,i, u_gf,uBarra_f, uBarra_g,theta);

                    U(n+1,i)= U(n,i) - 0.5*delta_t*(1.0/deltaX)*(uBarra_f* u_gf(1) - uBarra_g * u_gf(0)) + (delta_t*visc/ (deltaX*deltaX))*(U(n,i+1) - 2*U(n,i) + U(n,i-1));

                }
            }
            cout << "deu certo adb" <<endl;
            break;
        }
    }
    return U;
}

MatrixXf Burgers_Viscosa::InsereCondicoesIniciais(MatrixXf M, VectorXf X)
{
    for(int i=0; i<M.rows();i++)
    {
        for(int j=0;j<M.cols();j++)
        {
            if(i==0)
            {
                //APLICA CONDIÇÃO INICIAL
                M(i,j) = sin(2.0*pi*X(j));

            }else if(j==0)
            {
                //APLICA CONDIÇOES DE CONTORNO
                M(i,j)= this->contornoInicio;
            }else if(j==M.cols()-1)
            {
                //APLICA CONDIÇOES DE CONTORNO
                M(i,j)= this->contornoFim;
            }else
            {
                M(i,j)=0;
            }
        }
    }
    //ImprimeMatriz(M);
    return M;
}
