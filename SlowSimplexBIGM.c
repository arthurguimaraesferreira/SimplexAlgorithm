// Trabalho Prático 2
// Pesquisa Operacional (DCC035) - TZ 2025/1
// DCC/UFMG
// Aluno: Arthur Guimarães Ferreira
// Matrícula: 2023034781
// Professor: Marcio Costa Santos

// Algoritmo Simplex para resolver problemas do tipo:
// max c^t x
// st. Ax = b
//     x >= 0

// Método BigM de minimização:
// para resolver o problema de maximização,
// com esse algoritmo, fazemos a transformação abaixo.

// max c^t x        →       min (-c)^t x
// st. Ax = b               st. Ax = b
//     x >= 0                   x >= 0

// Resultado final = oposto do resultado!

// Videoaula explicativa do algoritmo: 
// https://www.youtube.com/watch?v=G6X2kzZvgJI&t=967s

#include <stdio.h>
double BIGM = -100000.0;

int main() {
    int ret = 0;
    double ret2 = 0.0;

    int n = 0;
    int m = 0;
    ret = scanf("%d %d", &n, &m);

    int l = n + 1;
    int c = n + m + 1;

    double tableau[l][c];
    double C[c-1];
    double Z[c-1];
    double custos_reduzidos[c-1];


    // Configuração F. Objetivo 
    for (int j = 0; j < m; j++) {
        ret2 = scanf("%lf", &tableau[0][j]);
        C[j] = tableau[0][j];
    }
    for (int k = m; k < c-1; k++) {
        tableau[0][k] = BIGM;
        C[k] = tableau[0][k];
    }
    tableau[0][c-1] = 0.0;

    // Leiturado Tableau
    for (int i = 1; i < l; i++) {
        for (int j = 0; j < m+1; j++) {
            if (j != m) {
                ret2 = scanf("%lf", &tableau[i][j]);
            } else {
                ret2 = scanf("%lf", &tableau[i][c-1]);
            }
        }

        // Configuando as variáveis auxiliares
        for (int k = m; k < c-1; k++) {
            if ((k-m) == (i-1)) {
                tableau[i][k] = 1.0;
            } else {
                tableau[i][k] = 0.0;
            }
        }
    }

    // Simplex com BIG-M

    int variaveis_base[n];
    for(int k = m; k < c-1; k++) {
        variaveis_base[k-m] = k;
    }


    while(1) {
        // 0. Calcular Custos Reduzidos
        for(int k = 0; k < c-1; k++) {
            custos_reduzidos[k] = 0;
            Z[k] = 0;
        }
        
        for (int k = 0; k < c-1; k++) {
            for (int i = 1; i < l; i++) {
                Z[k] = Z[k] + (C[variaveis_base[i-1]] * tableau[i][k]);
            }
            custos_reduzidos[k] = C[k] - Z[k];
        }

        // 1. Definir Coluna Pivô
        int col_pivo = -1;
        double maior_custo_reduzido = 0.0;
        for (int j = 0; j < c-1; j++) {
            if (custos_reduzidos[j] > maior_custo_reduzido) {
                maior_custo_reduzido = custos_reduzidos[j];
                col_pivo = j;
            }
        }
        if (col_pivo == -1) break; // Ótimo

        // 2. Definir Linha Pivô
        int lin_pivo = -1;
        double menor_divisao = 1e20;
        for (int i = 1; i < l; i++) {
            double aij = tableau[i][col_pivo];
            double bi  = tableau[i][c-1];
            
            if (aij > 0 && bi >= 0) {
                double razao = bi / aij;
                if (razao < menor_divisao) {
                    menor_divisao = razao;
                    lin_pivo = i;
                }
            }
        }
        if (lin_pivo == -1) {
            printf("ilimitado\n");
            return 0;
        }

        // 3. Mudar base
        variaveis_base[lin_pivo-1] = col_pivo;


        // 4. Pivoteamento
        double piv = tableau[lin_pivo][col_pivo];

        for (int k = 0; k < c; k++) {
            tableau[lin_pivo][k] = tableau[lin_pivo][k] / piv;
        }

        for (int i = 0; i < l; i++) {
            if (i != lin_pivo) {
                double f = tableau[i][col_pivo];
                for (int j = 0; j < c; j++) {
                    tableau[i][j] = tableau[i][j] - f*tableau[lin_pivo][j];
                }
            }
        }
    }

    int inviable = 0;
    for (int i = 0; i < n; i++) {
        int var = variaveis_base[i];
        if (var >= m && var < c-1) {
            double valor = tableau[i+1][c-1];
            if (valor > 1e-8) {        // tolerância para zeros numéricos
                inviable = 1;
                break;
            }
        }
    }
    if (inviable) {
        printf("inviavel\n");
        return 0;
    }

    printf("otimo\n%.3lf\n", -tableau[0][c-1]);

    return 0;
}