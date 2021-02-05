#include "upwind.h"

using namespace Eigen;
using namespace std;


Upwind::Upwind(float theta)
{
    //construtor
    cout<< "1 = FSLS " <<endl;
    cout<< "2 = ADBQUICKEST"<<endl;
    cin >> tipoDeUpwind;

    switch(tipoDeUpwind)
    {
    case(0):
    {
        cout<<"Erro ao escolher o tipo de upwind"<<endl;
        break;
    }
    case(1):
    {
        //cout<<"Valor de Beta para FSLS: "<<endl;
        //cin >> this->beta;
        beta= 1.5;
        break;
    }
    case(2):
    {
        float num_a = 2-3*abs(theta)+theta*theta;
        float den_a = 7 - 6*theta- 3*abs(theta)+ 2*theta*theta;
        this->a = num_a/den_a;

        float num_b = -4+6*theta- 3*abs(theta)+ theta*theta;
        float den_b = -5+ 6*theta - 3*abs(theta) + 2*theta*theta;
        this->b = num_b/den_b;
    }
    }
}

Vector2f Upwind::ADBQUICKEST(MatrixXf M, int i, int j, Vector2f v,float u_f, float u_g, float theta) //v[0] é u_g e v[1] é u_f /i=cont no tempo e j=cont no espaço
{

    //calcula as propriedade convectada
    if(u_g>=0.0)
    {
        //calcula as propriedade convectada qdo a velocidade de transporte é positiva
        if(j == 1)
        {
            v(0)= 0.5*(M(i,j-1)+M(i,j));
        }else
        {
            float phi_chapeu = (M(i,j-1)-M(i,j-2))/(M(i,j)-M(i,j-2));
            if(phi_chapeu<0 || 1<phi_chapeu)
            {
                v(0)= M(i,j-1);
            }else if(phi_chapeu>=0 && phi_chapeu<a)
            {
                float A = (2-theta)*phi_chapeu;
                v(0) = M(i,j-2) + (M(i,j)-M(i,j-2))*A;
            }else if(phi_chapeu>=a && phi_chapeu<=b)
            {
                float A = phi_chapeu+0.5*(1-abs(theta))-(1/6)*(1-theta*theta)*(1-2*phi_chapeu);
                v(0) = M(i,j-2) + (M(i,j)-M(i,j-2))*A;
            }else if(phi_chapeu>b && phi_chapeu<=1)
            {
                float A = 1-theta+theta*phi_chapeu;
                v(0) = M(i,j-2) + (M(i,j)-M(i,j-2))*A;
            }
        }

    }
    if(u_g<0.0)
    {
        //calcula as propriedade convectada qdo a velocidade de transporte é positiva
        if(fabs(M(i,j-1)-M(i,j+1)) <= threshold)
        {
            v(0) = M(i,j);
        }else
        {
            float phi_chapeu = (M(i,j)-M(i,j+1))/(M(i,j-1)-M(i,j+1));
            if(phi_chapeu<0 || 1<phi_chapeu)
            {
                v(0)= M(i,j);
            }else if(phi_chapeu>=0 && phi_chapeu<a)
            {
                float A = (2-theta)*phi_chapeu;
                v(0) = M(i,j+1) + (M(i,j-1)-M(i,j+1))*A;
            }else if(phi_chapeu>=a && phi_chapeu<=b)
            {
                float A = phi_chapeu+0.5*(1-abs(theta))-(1/6)*(1-theta*theta)*(1-2*phi_chapeu);
                v(0) = M(i,j+1) + (M(i,j-1)-M(i,j+1))*A;
            }else if(phi_chapeu>b && phi_chapeu<=1)
            {
                float A = 1-theta+theta*phi_chapeu;
                v(0) = M(i,j+1) + (M(i,j-1)-M(i,j+1))*A;
            }
        }

    }
    if(u_f>=0.0)
    {
        //calcula as propriedade convectada qdo a velocidade de transporte é positiva
        if(fabs(M(i,j+1)-M(i,j-1)) <= threshold)
        {
            v(1)= M(i,j);
        }else
        {
            float phi_chapeu = (M(i,j)-M(i,j-1))/(M(i,j+1)-M(i,j-1));
            if(phi_chapeu<0 || 1<phi_chapeu)
            {
                v(1)= M(i,j);
            }else if(phi_chapeu>=0 && phi_chapeu<a)
            {
                float A = (2-theta)*phi_chapeu;
                v(1) = M(i,j-1) + (M(i,j+1)-M(i,j-1))*A;
            }else if(phi_chapeu>=a && phi_chapeu<=b)
            {
                float A = phi_chapeu+0.5*(1-abs(theta))-(1/6)*(1-theta*theta)*(1-2*phi_chapeu);
                v(1) = M(i,j-1) + (M(i,j+1)-M(i,j-1))*A;
            }else if(phi_chapeu>b && phi_chapeu<=1)
            {
                float A = 1-theta+theta*phi_chapeu;
                v(1) = M(i,j-1) + (M(i,j+1)-M(i,j-1))*A;
            }
        }

    }
    if(u_f<0.0)
    {
        //calcula as propriedade convectada qdo a velocidade de transporte é positiva
        if(j == M.cols()-2)
        {
            v(1)= 0.5*(M(i,j-1)+M(i,j+1));
        }else
        {
            float phi_chapeu = (M(i,j+1)-M(i,j+2))/(M(i,j)-M(i,j+2));
            if(phi_chapeu<0 || 1<phi_chapeu)
            {
                v(1)= M(i,j+1);
            }else if(phi_chapeu>=0 && phi_chapeu<a)
            {
                float A = (2-theta)*phi_chapeu;
                v(1) = M(i,j+2) + (M(i,j)-M(i,j+2))*A;
            }else if(phi_chapeu>=a && phi_chapeu<=b)
            {
                float A = phi_chapeu+0.5*(1-abs(theta))-(1/6)*(1-theta*theta)*(1-2*phi_chapeu);
                v(1) = M(i,j+2) + (M(i,j)-M(i,j+2))*A;
            }else if(phi_chapeu>b && phi_chapeu<=1)
            {
                float A = 1-theta+theta*phi_chapeu;
                v(1) = M(i,j+2) + (M(i,j)-M(i,j+2))*A;
            }
        }

    }
    return v;
}

Vector2f Upwind::FSLS(MatrixXf M, int i, int j, Vector2f v,float u_f, float u_g)
{

    //verifica em qual caso a variavel convectada se encontra
    if(u_g>=0.0)
    {
        //calcula as propriedade convectada qdo a velocidade de transporte é positiva
        if(j == 1.0)
        {
            v(0)= 0.5*(M(i,j-1)+M(i,j));
        }else
        {
            if(fabs(M(i,j)-M(i,j-2))<=threshold)
            {
                v(0)= M(i,j-1);
            }else
            {
                float phi_chapeu = (M(i,j-1)-M(i,j-2))/(M(i,j)-M(i,j-2));
                if(phi_chapeu<0.0 || 1.0<phi_chapeu)
                {
                    v(0)= M(i,j-1);
                }else
                {
                    float A = (-2.0*beta+4.0)*pow(phi_chapeu,4.0) + (4.0*beta-8.0)*pow(phi_chapeu,3.0)+ ((-5.0*beta+8.0)/2.0)*pow(phi_chapeu,2.0)+ ((beta+2.0)/2.0)*phi_chapeu;
                    v(0) = M(i,j-2) + (M(i,j)-M(i,j-2))*A;
                }
            }

        }

    }else if(u_g<0.0)
    {
        //calcula as propriedade convectada qdo a velocidade de transporte é positiva
        if(fabs(M(i,j-1)-M(i,j+1)) <= threshold)
        {
            v(0)= M(i,j);
        }else
        {
            float phi_chapeu = (M(i,j)-M(i,j+1))/(M(i,j-1)-M(i,j+1));
            if(phi_chapeu<0.0 || 1.0<phi_chapeu)
            {
                v(0)= M(i,j);
            }else
            {
                float A = (-2.0*beta+4.0)*pow(phi_chapeu,4.0) + (4.0*beta-8.0)*pow(phi_chapeu,3.0)+ ((-5.0*beta+8.0)/2.0)*pow(phi_chapeu,2.0)+ ((beta+2.0)/2.0)*phi_chapeu;
                v(0) = M(i,j+1) + (M(i,j-1)-M(i,j+1))*A;
            }
        }
    }
    if(u_f>=0.0)
    {
        //calcula as propriedade convectada qdo a velocidade de transporte é positiva
        if(fabs(M(i,j+1)-M(i,j-1))<=threshold)
        {
            v(1)= M(i,j);
        }else
        {
            float phi_chapeu = (M(i,j)-M(i,j-1))/(M(i,j+1)-M(i,j-1));
            if(phi_chapeu<0 || 1<phi_chapeu)
            {
                v(1)= M(i,j);
            }else
            {
                float A = (-2.0*beta+4.0)*pow(phi_chapeu,4.0) + (4.0*beta-8.0)*pow(phi_chapeu,3.0)+ ((-5.0*beta+8.0)/2.0)*pow(phi_chapeu,2.0)+ ((beta+2.0)/2.0)*phi_chapeu;
                v(1) = M(i,j-1) + (M(i,j+1)-M(i,j-1))*A;
            }
        }
    }else if(u_f<0.0)
    {
        //calcula as propriedade convectada qdo a velocidade de transporte é positiva
        if(j == M.cols()-2)
        {
            v(1)= 0.5*(M(i,j)+M(i,j+1));
        }else
        {
            if(fabs(M(i,j)-M(i,j+2)) <= threshold)
            {
                v(1)= M(i,j+1);
            }else
            {
                float phi_chapeu = (M(i,j+1)-M(i,j+2))/(M(i,j)-M(i,j+2));
                if(phi_chapeu<0.0 || 1.0<phi_chapeu)
                {
                    v(1)= M(i,j+1);
                }else
                {
                    float A = (-2.0*beta+4.0)*pow(phi_chapeu,4.0) + (4.0*beta-8.0)*pow(phi_chapeu,3.0)+ ((-5.0*beta+8.0)/2.0)*pow(phi_chapeu,2.0)+ ((beta+2.0)/2.0)*phi_chapeu;
                    v(1) = M(i,j+2) + (M(i,j)-M(i,j+2))*A;
                }
            }
        }
    }
    return v;
}

