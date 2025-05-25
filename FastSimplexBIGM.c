// Trabalho Prático 2
// Pesquisa Operacional (DCC035) - TZ 2025/1
// DCC/UFMG
// Aluno: Arthur Guimarães Ferreira
// Professor: Marcio Costa Santos

// Algoritmo Simplex para resolver problemas do tipo:
// max c^t x
// st. Ax = b
//     x >= 0

// Método Big-M de *minimização*:
// para resolver o problema de maximização
// com esse algoritmo, fazemos a transformação abaixo.

// max c^t x        →       min (-c)^t x
// st. Ax = b               st. Ax = b
//     x >= 0                   x >= 0

// Resultado final = oposto do resultado!

// Videoaula explicativa do algoritmo: 
// https://www.youtube.com/watch?v=G6X2kzZvgJI&t=967s

// Nesta implementação, o objetivo é fazer o Simplex
// que rode O MAIS RÁPIDO POSSÍVEL!
// Uma versão mais lenta, porém mais simples do
// código se encontra em: https://github.com/arthurguimaraesferreira/SimplexAlgorithm
// Esta versão com foco em baixo tempo de execução foi
// melhorada a partir da versão acima.

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

double BIGM = -100000.0;


int main() {
    int ret = 0;
    double ret2 = 0.0;

    int n = 0;
    int m = 0;
    ret = scanf("%d %d", &n, &m);

    int l = n + 1;
    int c = n + m + 1;

    double * restrict tableau;
    //double * restrict tableau = malloc(sizeof(double) * l * c);
    if (posix_memalign((void**)&tableau, 64, sizeof(double) * l * c) != 0) {
        fprintf(stderr, "Erro de alocação alinhada\n");
        return 1;
    }
    #define tableau(i,j) tableau[(i) * c + (j)]

    double C[c];
    double Z[c];

    // Configuração F. Objetivo 
    for (int j = 0; j < m; j++) {
        ret2 = scanf("%lf", &tableau(0, j));
        C[j] = tableau(0, j);
        Z[j] = 0.0;
    }
    for (int k = m; k < c-1; k++) {
        tableau(0, k) = BIGM;
        C[k] = tableau(0, k);
        Z[k] = 0.0;
    }
    tableau(0, c-1) = 0.0;
    C[c-1] = 0.0;
    Z[c-1] = 0.0;

    // Leiturado Tableau
    for (int i = 1; i < l; i++) {
        for (int j = 0; j < m; j++) {
            ret2 = scanf("%lf", &tableau(i, j));
        }
        ret2 = scanf("%lf", &tableau(i, c-1));

        // Configuando as variáveis auxiliares
        for (int k = m; k < c-1; k++) {
            if ((k-m) == (i-1)) {
                tableau(i, k) = 1.0;
            } else {
                tableau(i, k) = 0.0;
            }
        }
    }
    
/*************************************************************************************** */
/*************************************************************************************** */
/*************************************************************************************** */

    // SIMPLEX, MÉTODO BIGM

    // Inicia as váriaveis base
    int variaveis_base[n];
    for(int k = m; k < c-1; k++) {
        variaveis_base[k-m] = k;
    }

    // 0. Calcula o primeiro custo reduzido. Os próximos são calculados
    //    dentro do loop do Simplex, pela linha tableau[0][k], 0 <= k < c
    for (int k = 0; k < c; k++) {
        for (int i = 1; i < l; i++) {
            Z[k] = Z[k] + (C[variaveis_base[i-1]] * tableau(i, k));
        }
        tableau(0, k) = C[k] - Z[k];
    }

    // Loop do Simplex
    while(1) {

        // 1. Definir Coluna Pivô
        int col_pivo = -1;
        double maior_custo_reduzido = 0.0;
        double *row0 = &tableau(0, 0);
        int j = 0;
        // unroll de 4 em 4
        for (; j + 3 < c-1; j += 4) {
            if (row0[j]   > maior_custo_reduzido) { maior_custo_reduzido = row0[j];   col_pivo = j;   }
            if (row0[j+1] > maior_custo_reduzido) { maior_custo_reduzido = row0[j+1]; col_pivo = j+1; }
            if (row0[j+2] > maior_custo_reduzido) { maior_custo_reduzido = row0[j+2]; col_pivo = j+2; }
            if (row0[j+3] > maior_custo_reduzido) { maior_custo_reduzido = row0[j+3]; col_pivo = j+3; }
        }
        for (; j < c-1; j++) {
            if (row0[j] > maior_custo_reduzido) {
                maior_custo_reduzido = row0[j];
                col_pivo = j;
            }
        }
        if (col_pivo == -1) break; // Ótimo


        // 2. Definir Linha Pivô (ponteiro para cada linha)
        int lin_pivo = -1;
        double menor_divisao = 1e20;
        for (int i = 1; i < l; i++) {
            double *rowi = &tableau(i, 0);
            double aij = rowi[col_pivo];
            double bi = rowi[c-1];
            if (aij > 0.0 && bi >= 0.0) {
                double razao = bi / aij;
                if (razao < menor_divisao) {
                    menor_divisao = razao;
                    lin_pivo = i;
                }
            }
        }
        if (lin_pivo == -1) {
            printf("ilimitada\n");
            free(tableau);
            return 0;
        }


        // 3. Mudar base
        variaveis_base[lin_pivo-1] = col_pivo;


        // 4. Pivoteamento
        double piv = tableau(lin_pivo, col_pivo);

        double *pivo_row = &tableau(lin_pivo, 0);
        double inv_piv = 1.0f / piv;
        for (int k = 0; k < c; k++) {
            pivo_row[k] = pivo_row[k] * inv_piv;
        }

        // Elimina a coluna pivô em duas fases, sem branch em i
        // e usando ponteiros de linha para remover cálculo de índices

        // fase 1: i = 0 … lin_pivo-1
        for (int i = 0; i < lin_pivo; i++) {
            double *row_i = &tableau(i, 0);
            double f = row_i[col_pivo];
            int j = 0;
            // unroll de 4 em 4
            for (; j + 3 < c; j += 4) {
                row_i[j]   -= f * pivo_row[j];
                row_i[j+1] -= f * pivo_row[j+1];
                row_i[j+2] -= f * pivo_row[j+2];
                row_i[j+3] -= f * pivo_row[j+3];
            }
            for (; j < c; j++) {
                row_i[j] -= f * pivo_row[j];
            }
        }

        // fase 2: i = lin_pivo+1 … l-1
        for (int i = lin_pivo + 1; i < l; i++) {
            double *row_i = &tableau(i, 0);
            double f = row_i[col_pivo];
            int j = 0;
            // unroll de 4 em 4
            for (; j + 3 < c; j += 4) {
                row_i[j]   -= f * pivo_row[j];
                row_i[j+1] -= f * pivo_row[j+1];
                row_i[j+2] -= f * pivo_row[j+2];
                row_i[j+3] -= f * pivo_row[j+3];
            }
            for (; j < c; j++) {
                row_i[j] -= f * pivo_row[j];
            }
        }


    }

    // Verificação de viabilidade ao final da execução do Simplex
    bool inviavel = false;
    for (int i = 0; i < n; i++) {
        if (variaveis_base[i] >= m && variaveis_base[i] < c-1) {
            double valor = tableau(i+1, c-1);
            if (valor > 1e-8) {        // tolerância para zeros numéricos
                inviavel = true;
                break;
            }
        }
    }
    if (inviavel) {
        printf("inviavel\n");
        free(tableau);
        return 0;
    }

    printf("otima\n%.3f\n", -tableau(0, c-1));

    free(tableau);
    return 0;
}