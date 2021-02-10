//#include <QCoreApplication>
#include "../Eigen/Eigen/Dense"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "burgers_viscosa.h"
#include "adveccao_linear.h"
#include "upwind.h"


using namespace Eigen;
using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------------------------*/
void ImprimeMatriz(MatrixXf M)
{
    for(int i=0; i<M.rows();i++)
    {
        for(int j=0;j<M.cols();j++)
        {
            cout << M(i,j) <<" ";
        }
        cout << endl;
    }
}

void SalvaArquivo(MatrixXf M, FILE *arquivo, VectorXf tempo, VectorXf X)
{
    arquivo = fopen("MatrizSolucaoTrabalho1.dat","w");
    if(!arquivo)
    {
        printf("\nerro ao abrir arquivo");
    }else
    {
        for(int i=0;i<M.rows();i++)
        {
            for(int j=0;j<M.cols();j++)
            {
               fprintf(arquivo, "%f  %f  %f  \n", tempo(i), X(j), M(i,j));
            }
        }
        cout << "arquivo atualizado" <<endl;
        fclose(arquivo);
    }
}

void SalvaArquivoSemTempo(MatrixXf M, FILE *arquivo, VectorXf X)
{
    arquivo = fopen("SolucaoTrabalho1SemTempo.dat","w");
    if(!arquivo)
    {
        printf("\nerro ao abrir arquivo");
    }else
    {
        int i = M.rows()-1;
        for(int j=0;j<M.cols();j++)
        {
           fprintf(arquivo, "%f  \n", M(i,j));
        }
        cout << "arquivo atualizado" <<endl;
        fclose(arquivo);
    }
}


//int main(int argc, char *argv[])
int main()
{
    //QCoreApplication a(argc, argv);

    FILE *arquivo;


    int qualEquacao = 0;
    while(true)
    {
        cout<<"Qual equacao deseja resolver?"<<endl;
        cout<<"1 - Burgers Viscosa"<<endl;
        cout<<"2 - Adveccao Linear"<<endl;
        cin >> qualEquacao;

        if(qualEquacao == 1)
        {
            Burgers_Viscosa Equacao;

            Equacao.contornoInicio = 0;
            Equacao.contornoFim = 0;

            Equacao.cont_t= floor(Equacao.tempo/Equacao.delta_t)+1;
            Equacao.contX = floor(Equacao.comprimento/Equacao.deltaX)+1;



            VectorXf X(Equacao.contX);
            VectorXf tempo(Equacao.cont_t);

            for(int c=0;c<X.rows();c++)
            {
                //inicia o vetor X com a posição de cada nó do dominio 1D com a borda esquerda sendo a coordenada 0
                if(c==0)
                {
                    X[c]=0;

                }else
                {
                    X[c]=X[c-1]+Equacao.deltaX;

                }
            }

            for(int c=0;c<tempo.rows();c++)
            {
                //inicia o vetor tempo com a posição de cada nó no decorrer do tempo
                if(c==0)
                {
                    tempo[c]=0;

                }else
                {
                    tempo[c]=tempo[c-1]+Equacao.delta_t;

                }
            }

            MatrixXf solucao = Equacao.CalculaEquacao(X);
            cout <<solucao.rows()<<"x"<<solucao.cols()<<endl;
            //ImprimeMatriz(solucao);

            SalvaArquivo(solucao, arquivo, tempo, X);
            SalvaArquivoSemTempo(solucao, arquivo, X);
        }
        if(qualEquacao == 2)
        {
            Adveccao_Linear Equacao;

            Equacao.contornoInicio = 0;
            Equacao.contornoFim = 0;

            Equacao.cont_t= floor(Equacao.tempo/Equacao.delta_t)+1;
            Equacao.contX = floor(Equacao.comprimento/Equacao.deltaX)+1;



            VectorXf X(Equacao.contX);
            VectorXf tempo(Equacao.cont_t);

            for(int c=0;c<X.rows();c++)
            {
                //inicia o vetor X com a posição de cada nó do dominio 1D com a borda esquerda sendo a coordenada 0
                if(c==0)
                {
                    X[c]=0;

                }else
                {
                    X[c]=X[c-1]+Equacao.deltaX;

                }
            }

            for(int c=0;c<tempo.rows();c++)
            {
                //inicia o vetor tempo com a posição de cada nó no decorrer do tempo
                if(c==0)
                {
                    tempo[c]=0;

                }else
                {
                    tempo[c]=tempo[c-1]+Equacao.delta_t;

                }
            }

            MatrixXf solucao = Equacao.CalculaEquacao(X);
            cout <<solucao.rows()<<"x"<<solucao.cols()<<endl;
            //ImprimeMatriz(solucao);

            SalvaArquivo(solucao, arquivo, tempo, X);
            SalvaArquivoSemTempo(solucao, arquivo, X);
        }
    }
    //return a.exec();
    return 0;
}
