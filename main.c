/*
Trabalho de Ajuste de Curvas - T2
Data: 27/09/2024
Integrantes: Bernardo Krause Rodrigues, Caio Freire Silva, Gabriel Tetzner Menegueti
*/

// Bibliotecas necessarias
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Tamanho dos dados
#define N 48

// Funcao que copia uma matriz para uma outra
void copiaMatriz(double **origem, double **destino, int tamL, int tamC) {
  for (int i = 0; i < tamL; i++) {
    for (int j = 0; j < tamC; j++) {
      destino[i][j] = origem[i][j];
    }
  }
}

// Funcao que copia um vetor para um outro
void copiaVetor(double *origem, double *destino, int tam) {
  for (int i = 0; i < tam; i++) {
    destino[i] = origem[i];
  }
}

// Retorna o valor do Multiplicador
double calculaM(double a, double b) { return a / b; }

// Realiza a substituicao sucessiva
void realizaSubstituicaoSucessiva(double **A, double *X, double *B,
                                  int dimensao) {
  int i, j;
  for (i = dimensao - 1; i >= 0; i--) {
    for (j = i + 1; j < dimensao; j++) {
      B[i] -= A[i][j] * X[j];
    }
    X[i] = B[i] / A[i][i];
  }
}

// Realiza as etapas da eliminacao de gauss
void realizaEtapa(double **A, double *B, int dimensao) {
  int k, l, c;
  double multiplicador;
  for (k = 0; k < dimensao - 1; k++) {
    for (l = k + 1; l < dimensao; l++) {
      multiplicador = calculaM(A[l][k], A[k][k]);
      for (c = 0; c < dimensao; c++) {
        A[l][c] -= multiplicador * A[k][c];
      }
      B[l] -= multiplicador * B[k];
    }
  }
}

// Eliminacao de Gauss
void eliminacaoGauss(double **A, double *X, double *B, int dimensao) {
  realizaEtapa(A, B, dimensao);
  realizaSubstituicaoSucessiva(A, X, B, dimensao);
}

// Funcao para realizar os somatorios necessarios para as devidas regressoes
void realizaSomatorios(double *X, double *Y, double *somatorio_1,
                       double *somatorio_x, double *somatorio_y,
                       double *somatorio_y2, double *somatorio_x2,
                       double *somatorio_x3, double *somatorio_x4,
                       double *somatorio_xy, double *somatorio_x2y,
                       double *somatorio_ylog, double *somatorio_xylog) {

  for (int i = 0; i < N; i++) {

    *somatorio_1 += 1.0;

    *somatorio_x += X[i];

    *somatorio_y += Y[i];

    *somatorio_y2 += Y[i] * Y[i];

    *somatorio_x2 += (X[i] * X[i]);

    *somatorio_x3 += (X[i] * X[i] * X[i]);

    *somatorio_x4 += (X[i] * X[i] * X[i] * X[i]);

    *somatorio_xy += (X[i] * Y[i]);

    *somatorio_x2y += (X[i] * X[i]) * (Y[i]);

    *somatorio_ylog += log(Y[i]);

    *somatorio_xylog += (X[i] * log(Y[i]));
  }
}

// Funcao para calcular o coeficiente de determinacao de cada regressao
double calculaCoeficienteDeterminacao(double somatorio_yu2, double somatorio_y2,
                                      double somatorio_y, int qtdPontos) {
  return 1 - (somatorio_yu2 /
              (somatorio_y2 - ((somatorio_y * somatorio_y) / qtdPontos)));
}

// Imprime resposta ira imprimir os coeficientes de determinacao e qual curva se ajusta melhor
void imprimeResposta(double r2_linear, double r2_quadratico,
                     double r2_exponencial) {
  printf("\nCoeficientes de determinacao: \n\nLinear: %lf\n\nQuadratico: "
         "%lf\n\nExponencial: %lf\n\n",
         r2_linear, r2_quadratico, r2_exponencial);

  if (r2_linear > r2_quadratico && r2_linear > r2_exponencial) {
    printf("Melhor ajuste: Regressao Linear");
  } else if (r2_quadratico > r2_linear && r2_quadratico > r2_exponencial) {
    printf("Melhor ajuste: Regressao Quadratica");
  } else {
    printf("Melhor ajuste: Regressao Exponencial");
  }
}

// Main
int main(void) {

  int i;

  // Declaracao do vetor X (anos t)
  double X[N] = {77,   119,  205,  260,  343,  415,  425,  438,  502,  580,
                  604,  675,  696,  770,  802,  822,  897,  965,  970,  1027,
                  1094, 1156, 1192, 1282, 1345, 1405, 1429, 1493, 1516, 1597,
                  1678, 1721, 1724, 1812, 1873, 1947, 1950, 2012, 2047, 2127,
                  2153, 2157, 2210, 2298, 2332, 2358, 2449, 2503};

  // Declaracao do vetor Y (quantidade N)
  double Y[N] = {
      50870643080, 46297918240, 38282822421, 34080460561, 28088573347,
      24175635810, 23588757299, 22718540971, 19736321394, 16655820340,
      15956176275, 13458062313, 13121778454, 10939070444, 10054385447,
      9824068939,  8451625483,  7251508116,  7118718866,  6316609405,
      5464716073,  4511225310,  4322167723,  3709676790,  3226060519,
      2761872970,  2553745475,  2220343372,  2046714895,  1667362442,
      1648561844,  1167013186,  1211599434,  1291077808,  1099297084,
      791212548,   662664385,   721592837,   501015203,   536033559,
      510953997,   583955494,   403003225,   574150933,   412508389,
      131844628,   330342316,   326167250};

  // Declaracao das matrizes, vetores e das variaveis para armazenar os somatorios
  double **A, **auxA, *betas, *B, *auxB,
      somatorio_1 = 0.0, somatorio_x = 0.0, somatorio_y = 0.0,
      somatorio_y2 = 0.0, somatorio_x2 = 0.0, somatorio_x3 = 0.0,
      somatorio_x4 = 0.0, somatorio_xy = 0.0, somatorio_x2y = 0.0,
      somatorio_yu2 = 0.0, somatorio_ylog = 0.0, somatorio_xylog = 0.0, r2_linear, r2_quadratico, r2_exponencial;

  // Alocacao dinamica das matrizes e vetores
  A = (double **)malloc(3 * sizeof(double *));
  auxA = (double **)malloc(3 * sizeof(double *));
  
  for (int i = 0; i < 3; i++) {
    (A)[i] = (double *)malloc(3 * sizeof(double));
    (auxA)[i] = (double *)malloc(3 * sizeof(double));
  }

  betas = (double *)malloc(3 * sizeof(double));
  
  B = (double *)malloc(3 * sizeof(double));
  auxB = (double *)malloc(3 * sizeof(double));

  // Passa as variaveis dos somatorios como referencia e realiza os somatorios necessarios para as regressoes
  realizaSomatorios(X, Y, &somatorio_1, &somatorio_x, &somatorio_y,
                    &somatorio_y2, &somatorio_x2, &somatorio_x3, &somatorio_x4,
                    &somatorio_xy, &somatorio_x2y, &somatorio_ylog,
                    &somatorio_xylog);

  // Montar sistema linear a ser resolvido

  A[0][0] = somatorio_1;
  A[0][1] = somatorio_x;
  A[0][2] = somatorio_x2;
  A[1][0] = somatorio_x;
  A[1][1] = somatorio_x2;
  A[1][2] = somatorio_x3;
  A[2][0] = somatorio_x2;
  A[2][1] = somatorio_x3;
  A[2][2] = somatorio_x4;
  B[0] = somatorio_y;
  B[1] = somatorio_xy;
  B[2] = somatorio_x2y;
  
  // Preserva a matriz e o vetor copiando para os seus respectivos auxiliares
  copiaMatriz(A, auxA, 3, 3);
  copiaVetor(B, auxB, 3);

  // Regressao Linear

  // Realiza a eliminacao de gauss para resolver o sistema linear
  eliminacaoGauss(auxA, betas, auxB, 2);

  // Imprime os coeficientes da regressao
  printf("\n\nCoeficientes da regressao linear:\n");
  printf("b0 = %lf\n", betas[0]);
  printf("b1 = %lf\n", betas[1]);

  // Calcula o somatorio de yu2 com a funcao linear calculada
  for (i = 0; i < N; i++)
    somatorio_yu2 += (Y[i] - (betas[0] + (betas[1] * X[i]))) *
                     (Y[i] - (betas[0] + (betas[1] * X[i])));
  
  // Calcula coeficiente de determinacao
  r2_linear = calculaCoeficienteDeterminacao(somatorio_yu2, somatorio_y2,
                                                    somatorio_y, N);
  // Regressao Quadratica

  // Preserva a matriz e o vetor copiando para os seus respectivos auxiliares
  copiaMatriz(A, auxA, 3, 3);
  copiaVetor(B, auxB, 3);

  // Realiza a eliminacao de gauss para resolver o sistema linear
  eliminacaoGauss(auxA, betas, auxB, 3);

  // Imprime os coeficientes da regressao
  printf("\n\nCoeficientes da regressao quadratica:\n");
  printf("b0 = %lf\n", betas[0]);
  printf("b1 = %lf\n", betas[1]);
  printf("b2 = %lf\n", betas[2]);

  // Reseta somatorio de yu2
  somatorio_yu2 = 0.0;

  // Calcula o somatorio de yu2 com a funcao quadratica calculada
  for (i = 0; i < N; i++)
    somatorio_yu2 +=
        (Y[i] - (betas[0] + (betas[1] * X[i]) + (betas[2] * (X[i] * X[i])))) *
        (Y[i] - (betas[0] + (betas[1] * X[i]) + (betas[2] * (X[i] * X[i]))));

  // Calcula coeficiente de determinacao
  r2_quadratico = calculaCoeficienteDeterminacao(
      somatorio_yu2, somatorio_y2, somatorio_y, N);

  // Regressao Exponencial

  // Preserva a matriz e o vetor copiando para os seus respectivos auxiliares
  copiaMatriz(A, auxA, 3, 3);
  copiaVetor(B, auxB, 3);

  // Coloca os somatorios de log y e xy no vetor B
  auxB[0] = somatorio_ylog;
  auxB[1] = somatorio_xylog;

  // Realiza a eliminacao de gauss para resolver o sistema linear
  eliminacaoGauss(auxA, betas, auxB, 2);

  // Imprime os coeficientes da regressao
  printf("\n\nCoeficientes da regressao exponencial:\n");
  printf("b0 = %lf\n", exp(betas[0]));
  printf("b1 = %lf\n", betas[1]);
  
  // Reseta somatorio de yu2
  somatorio_yu2 = 0.0;

  // Calcula o somatorio de yu2 com a funcao exponencial calculada
  
  for (i = 0; i < N; i++) somatorio_yu2 += ((Y[i]) - (exp(betas[0]) * exp(betas[1] * X[i]))) *
                     ((Y[i]) - (exp(betas[0]) * exp(betas[1] * X[i])));
  

  // Calcula coeficiente de determinacao
  r2_exponencial = calculaCoeficienteDeterminacao(
      somatorio_yu2, somatorio_y2, somatorio_y, 48);

  // Imprime os coeficientes de determinacao e qual curva se ajusta melhor
  imprimeResposta(r2_linear, r2_quadratico, r2_exponencial);

  return 0;
}