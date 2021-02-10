#include "adveccao_linear.h"
#include "upwind.h"
#include <stdlib.h>

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
    this->comprimento = 2;
    //cout<<"Passo em X: "<<endl;
    //cin >> this->deltaX;
    this->deltaX=0.005;
    cout<<"Tempo total: "<<endl;
    cin >> this->tempo;
    //this->tempo = 0.5
    //cout<<"Passo do tempo: "<<endl;
    //cin >> this->delta_t;
    this->delta_t=theta*deltaX;

}

MatrixXf Adveccao_Linear::InsereCondicoesIniciais(MatrixXf M, VectorXf X)
{
    for(int i=0; i<M.rows();i++)
    {
        for(int j=0;j<M.cols();j++)
        {
            if(i==0)
            {
                //APLICA CONDIÇÃO INICIAL
                if ((X(j) >= 0.0) && (X(j) < 0.2))
                 {
                    M(i,j) = exp(-log(50.0)*((X(j)-0.15)/(0.05))*((X(j)-0.15)/(0.05)));
                 }
                 else if ((X(j) > 0.3) && (X(j) < 0.4))
                 {
                    M(i,j) = 1.0;
                 }
                 else if ((X(j) > 0.5) && (X(j) < 0.55))
                 {
                    M(i,j) = 20.0*X(j)-10.0;
                 }
                 else if ((X(j) >= 0.55) && (X(j) <0.6))
                 {
                    M(i,j) = 12.0-20.0*X(j);
                 }
                 else if ((X(j) > 0.7) && (X(j) < 0.8))
                 {
                    M(i,j) = sqrt(1.0-((X(j)-0.75)/(0.05))*((X(j)-0.75)/(0.05)));
                }
                 else if ((X(j) > 0.8) && (X(j) <= 2.0))
                 {
                    M(i,j) = 0.0;
                 }

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

MatrixXf Adveccao_Linear::CalculaEquacao(VectorXf X)
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
            double t = 0;
            for(int n=0; n<U.rows()-1;n++)
            {
                t += delta_t;
                system("CLS");
                cout << t << endl;
                //cout << delta_t << endl;
                for(int i=1; i<U.cols()-1;i++)
                {
                    //loop do espaço
                    uBarra_f = 2;
                    uBarra_g = 2;

                    u_gf = upwind.FSLS(U,n,i, u_gf,uBarra_f,uBarra_g);
                    float p1 = 0.5*delta_t*(1.0/deltaX)*(uBarra_f* u_gf(1) - uBarra_g * u_gf(0));
                    //float p2 = (delta_t*visc/ (deltaX*deltaX))*(U(n,i+1)-2.0*U(n,i) + U(n,i-1));
                    U(n+1,i)= U(n,i) - p1;

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
                    uBarra_f = 2;
                    uBarra_g = 2;

                    u_gf = upwind.ADBQUICKEST(U,n,i, u_gf,uBarra_f, uBarra_g,theta);

                    U(n+1,i)= U(n,i) - 0.5*delta_t*(1.0/deltaX)*(uBarra_f* u_gf(1) - uBarra_g * u_gf(0)); //+ (delta_t*visc/ (deltaX*deltaX))*(U(n,i+1) - 2*U(n,i) + U(n,i-1));

                }
            }
            cout << "deu certo adb" <<endl;
            break;
        }
    }
    return U;
}
