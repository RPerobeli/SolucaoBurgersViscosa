clear all;
close all;
clc;

tam = 101;
% entrada = fopen('SolucaoTrabalho1SemTempo.dat','r');
entrada = fopen('G:\Meu Drive\MESTRADO\Trimestre 3\TopMecanicaComputacional\Parte 2 - Hemodinâmica\SolucaoBurgersViscosa\build-TrabalhoUpwind-Desktop_Qt_5_11_1_MinGW_32bit-Debug\SolucaoTrabalho1SemTempo.dat','r');
% entrada = fopen('G:\Meu Drive\MESTRADO\Trimestre 3\TopMecanicaComputacional\Parte 2 - Hemodinâmica\SolucaoBurgersViscosa\build-TrabalhoUpwind-Desktop_Qt_5_11_1_MinGW_32bit-Debug\MatrizSolucaoTrabalho1.dat','r');
Resultado = fscanf(entrada,'%f',[tam,1]);

xi = 0;
xf = 1;
nel = length(Resultado)-1;
h = (xf- xi)/nel;

x = xi:h:xf;

plot(x,Resultado);