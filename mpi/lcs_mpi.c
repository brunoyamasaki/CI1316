#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

#define TAMANHO_ALFABETO 256

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
    int rank, num_procs;
    char *seqA = NULL, *seqB = NULL;
    int sizeA, sizeB;
    int **P_Matrix = NULL;
    double start_time, stop_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (rank == 0) {
        if (argc != 3) {
            printf("Uso: mpirun -np <N> %s <arquivo_sequencia_A> <arquivo_sequencia_B>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("Executando com %d processos.\n", num_procs);

        seqA = read_seq(argv[1], &sizeA);
        seqB = read_seq(argv[2], &sizeB);

        P_Matrix = (int **)malloc(TAMANHO_ALFABETO * sizeof(int *));
        for (int k = 0; k < TAMANHO_ALFABETO; k++) {
            P_Matrix[k] = (int *)calloc((sizeB + 1), sizeof(int));
        }
        calc_P_matrix(P_Matrix, seqB, sizeB);
    }

    start_time = MPI_Wtime();

    MPI_Bcast(&sizeA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sizeB, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        seqA = (char *)malloc((sizeA + 1) * sizeof(char));
        seqB = (char *)malloc((sizeB + 1) * sizeof(char));
        P_Matrix = (int **)malloc(TAMANHO_ALFABETO * sizeof(int *));
        for (int k = 0; k < TAMANHO_ALFABETO; k++) {
            P_Matrix[k] = (int *)calloc((sizeB + 1), sizeof(int));
        }
    }

    MPI_Bcast(seqA, sizeA + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(seqB, sizeB + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    for (int k = 0; k < TAMANHO_ALFABETO; k++) {
        MPI_Bcast(P_Matrix[k], sizeB + 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    int *linha_anterior = (int *)calloc(sizeB + 1, sizeof(int));
    int *linha_atual_global = (int *)calloc(sizeB + 1, sizeof(int));


    //aloca os vetores para o Allgatherv
    int *recv_counts = malloc(num_procs * sizeof(int));
    int *displs = malloc(num_procs * sizeof(int));

    //Calcula os chunks de trabalho para as colunas
    int items_per_proc = sizeB / num_procs;
    int remainder = sizeB % num_procs;
    int current_displ = 0;

    for (int p = 0; p < num_procs; p++) {
        recv_counts[p] = items_per_proc;
        if (p < remainder) {
            recv_counts[p]++;
        }
        displs[p] = current_displ;
        current_displ += recv_counts[p];
    }

    //Determina o chunk local deste processo
    int my_chunk_size = recv_counts[rank];
    int my_start_col = displs[rank] + 1; 
    int my_end_col = my_start_col + my_chunk_size;
    int* linha_atual_local = (int *)calloc(my_chunk_size, sizeof(int));




    for (int i = 1; i < sizeA + 1; i++) {
        char char_a = seqA[i - 1];
        int char_idx = (int)char_a;

        //Cada processo calcula seu chunk
        for (int j = my_start_col; j < my_end_col; j++) {
            int local_j = j - my_start_col; 
            if (char_a == seqB[j - 1]) {
                linha_atual_local[local_j] = linha_anterior[j - 1] + 1;
            } else {
                int p_val = P_Matrix[char_idx][j];
                if (p_val == 0) {
                    linha_atual_local[local_j] = linha_anterior[j];
                } else {
                    linha_atual_local[local_j] = max(linha_anterior[j], linha_anterior[p_val - 1] + 1);
                }
            }
        }

        MPI_Allgatherv(linha_atual_local,           
                      my_chunk_size,               
                      MPI_INT,
                      linha_atual_global + 1,      
                      recv_counts,                 
                      displs,                      
                      MPI_INT,
                      MPI_COMM_WORLD);

        linha_atual_global[0] = 0; 

        memcpy(linha_anterior, linha_atual_global, (sizeB + 1) * sizeof(int));
    }

    stop_time = MPI_Wtime();

    if (rank == 0) {
        int score = linha_atual_global[sizeB];
        printf("Tamanho da LCS: %d\n", score);
        printf("Tempo total de execucao (paralelo): %f segundos\n", stop_time - start_time);
    }


    free(seqA);
    free(seqB);
    for (int k = 0; k < TAMANHO_ALFABETO; k++) {
        free(P_Matrix[k]);
    }
    free(P_Matrix);
    free(linha_anterior);
    free(linha_atual_global);
    free(linha_atual_local);
    free(recv_counts); 
    free(displs);      

    MPI_Finalize();
    return 0;
}