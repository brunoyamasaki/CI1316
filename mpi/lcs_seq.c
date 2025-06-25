#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // Adicionado para medição de tempo

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

#define TAMANHO_ALFABETO 256

// Função para ler uma sequência de um arquivo
char* read_seq(char *fname, int* tamanho) {
    FILE *fseq = NULL;
    long size = 0;
    char *seq = NULL;
    int i = 0;
    fseq = fopen(fname, "rt");
    if (fseq == NULL) { printf("Erro lendo o arquivo %s\n", fname); exit(1); }
    fseek(fseq, 0L, SEEK_END);
    size = ftell(fseq);
    rewind(fseq);
    seq = (char *) calloc(size + 1, sizeof(char));
    if (seq == NULL) { printf("Erro alocando memoria para a sequencia %s.\n", fname); exit(1); }
    while (!feof(fseq)) {
        char c = fgetc(fseq);
        if ((c != '\n') && (c != EOF)) {
            seq[i] = c;
            i++;
        }
    }
    seq[i] = '\0';
    *tamanho = i;
    fclose(fseq);
    return seq;
}

// Função para calcular a Matriz P
void calc_P_matrix(int **P, char *seqB, int sizeB) {
    for (int i = 0; i < TAMANHO_ALFABETO; i++) {
        for (int j = 1; j < sizeB + 1; j++) {
            if (seqB[j - 1] == (char)i) {
                P[i][j] = j;
            } else {
                P[i][j] = P[i][j - 1];
            }
        }
    }
}


int main(int argc, char **argv) {
    char *seqA = NULL, *seqB = NULL;
    int sizeA, sizeB;
    int **P_Matrix = NULL;
    clock_t start_time, stop_time;

    if (argc != 3) {
        printf("Uso: %s <arquivo_sequencia_A> <arquivo_sequencia_B>\n", argv[0]);
        return 1;
    }

    printf("Executando de forma sequencial.\n");

    // Lógica que antes era do processo mestre, agora é a única lógica
    seqA = read_seq(argv[1], &sizeA);
    seqB = read_seq(argv[2], &sizeB);

    P_Matrix = (int **)malloc(TAMANHO_ALFABETO * sizeof(int *));
    for (int k = 0; k < TAMANHO_ALFABETO; k++) {
        P_Matrix[k] = (int *)calloc((sizeB + 1), sizeof(int));
    }
    
    start_time = clock(); // Inicia a medição de tempo

    calc_P_matrix(P_Matrix, seqB, sizeB);

    // Aloca os vetores para a linha anterior e a linha atual
    int *linha_anterior = (int *)calloc(sizeB + 1, sizeof(int));
    int *linha_atual = (int *)calloc(sizeB + 1, sizeof(int));

    // Loop principal sobre as linhas (caracteres de seqA)
    for (int i = 1; i < sizeA + 1; i++) {
        char char_a = seqA[i - 1];
        int char_idx = (int)char_a;

        // Lógica de cálculo sequencial da linha inteira
        for (int j = 1; j < sizeB + 1; j++) {
            if (char_a == seqB[j - 1]) {
                linha_atual[j] = linha_anterior[j - 1] + 1;
            } else {
                int p_val = P_Matrix[char_idx][j];
                int val1 = linha_anterior[j];
                int val2 = (p_val == 0) ? 0 : linha_anterior[p_val - 1] + 1;
                linha_atual[j] = max(val1, val2);
            }
        }
        
        // Prepara para a próxima iteração: a linha atual se torna a anterior
        // A chamada MPI_Allgatherv foi removida e este memcpy é tudo que precisamos.
        memcpy(linha_anterior, linha_atual, (sizeB + 1) * sizeof(int));
    }
    
    stop_time = clock(); // Finaliza a medição de tempo

    int score = linha_atual[sizeB];
    printf("Tamanho da LCS: %d\n", score);
    
    double tempo_total = ((double)(stop_time - start_time)) / CLOCKS_PER_SEC;
    printf("Tempo total de execucao (sequencial): %f segundos\n", tempo_total);
    
    // Libera toda a memória alocada
    free(seqA);
    free(seqB);
    for (int k = 0; k < TAMANHO_ALFABETO; k++) {
        free(P_Matrix[k]);
    }
    free(P_Matrix);
    free(linha_anterior);
    free(linha_atual);

    return 0;
}