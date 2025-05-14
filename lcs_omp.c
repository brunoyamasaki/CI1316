#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char* read_seq(char *fname) {
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	long size = 0;
	//sequence pointer
	char *seq = NULL;
	//sequence index
	int i = 0;

	//open file
	fseq = fopen(fname, "rt");
	if (fseq == NULL ) {
		printf("Error reading file %s\n", fname);
		exit(1);
	}

	//find out sequence size to allocate memory afterwards
	fseek(fseq, 0L, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);

	//allocate memory (sequence)
	seq = (char *) calloc(size + 1, sizeof(char));
	if (seq == NULL ) {
		printf("Erro allocating memory for sequence %s.\n", fname);
		exit(1);
	}

	//read sequence from file
	while (!feof(fseq)) {
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF))
			i++;
	}
	//insert string terminator
	seq[i] = '\0';

	//close file
	fclose(fseq);

	//return sequence pointer
	return seq;
}

mtype ** allocateScoreMatrix(int sizeA, int sizeB) {
	int i;
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) malloc((sizeB + 1) * sizeof(mtype *));
	for (i = 0; i < (sizeB + 1); i++)
		scoreMatrix[i] = (mtype *) malloc((sizeA + 1) * sizeof(mtype));
	return scoreMatrix;
}

void initScoreMatrix(mtype ** scoreMatrix, int sizeA, int sizeB) {
	int i, j;
	//Fill first line of LCS score matrix with zeroes
	for (j = 0; j < (sizeA + 1); j++)
		scoreMatrix[0][j] = 0;

	//Do the same for the first collumn
	for (i = 1; i < (sizeB + 1); i++)
		scoreMatrix[i][0] = 0;
}


int idx_char(char *seq, char x, int size) {
	int i;
	for (i = 0; i < size; i++) {
		if (seq[i] == x) {
			return i;
		}
	}
	fprintf(stderr, "Aviso: caractere não encontrado: '%c' (%d)\n", x, (int)x);
	return -1;
}

void calcula_p(mtype ** p_matrix, char * seqB, int sizeB, char * seqC, int sizeC) {
	#pragma omp para for
	for(int i = 0; i < sizeC ; i++) {
		for (int j = 0; j < sizeB + 1; j++) {
				if(j == 0){
					p_matrix[i][j] = 0;
				} else if(seqB[j - 1] == seqC[i]){
					p_matrix[i][j] = j;
				} else {
					p_matrix[i][j] = p_matrix[i][j - 1];
				}
			}
	}
}


int LCS_par(mtype ** scoreMatrix, mtype **p, int sizeA, int sizeB, int sizeC,char *seqA, char *seqB, char *seqC) {
	int i, j, t, s, idx_c;
	int *linha_ant = (int *)malloc((sizeB+1) * sizeof(int));

	for (i = 1; i < sizeA + 1; i++) {
		//idx_c = idx_char(seqC, seqA[i], sizeC);
		//if (idx_c == -1) {
		//	fprintf(stderr, "Erro: caractere não encontrado: '%c' (%d)\n", seqA[i - 1], (int)seqA[i - 1]);
		//	exit(EXIT_FAILURE);
			//continue;
		//}
		//for(int m = 0; m < sizeC - 1; m++){
		//	if (seqC[m] == seqA[i])
		//		idx_c = m;
		//}
		#pragma omp parallel for private(t,s) schedule(static)
		for (j = 1; j < sizeB + 1; j++) {
			if (p[idx_c][j] > 0) {
				t = 1;
			} else {
				t = 0;
			}
			if(t == 1){
				if (scoreMatrix[i-1][j] == scoreMatrix[i-1][p[idx_c][j]]) {
					s = 1;
				} else {
					s = 0;
				}
			} else {
				s = 1;
			}
			if (t == 0 || s == 0){
				scoreMatrix[i][j] = scoreMatrix[i-1][j];
			} else {
				scoreMatrix[i][j] = scoreMatrix[i-1][p[idx_c][j] - 1] + 1;
			}
		}
	}
	return scoreMatrix[sizeB][sizeA];
}

int LCS(mtype ** scoreMatrix, int sizeA, int sizeB, char * seqA, char *seqB) {
	int i, j;
	for (i = 1; i < sizeB + 1; i++) {
		for (j = 1; j < sizeA + 1; j++) {
			if (seqA[j - 1] == seqB[i - 1]) {
				/* if elements in both sequences match,
				 the corresponding score will be the score from
				 previous elements + 1*/
				scoreMatrix[i][j] = scoreMatrix[i - 1][j - 1] + 1;
			} else {
				/* else, pick the maximum value (score) from left and upper elements*/
				scoreMatrix[i][j] = max(scoreMatrix[i-1][j], scoreMatrix[i][j-1]);
			}
		}
	}
	return scoreMatrix[sizeB][sizeA];
}
void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix, int sizeA,
		int sizeB) {
	int i, j;

	//print header
	printf("Score Matrix:\n");
	printf("========================================\n");

	//print LCS score matrix allong with sequences

	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeA; j++)
		printf("%5c   ", seqA[j]);
	printf("\n");
	for (i = 0; i < sizeB + 1; i++) {
		if (i == 0)
			printf("    ");
		else
			printf("%c   ", seqB[i - 1]);
		for (j = 0; j < sizeA + 1; j++) {
			printf("%5d   ", scoreMatrix[i][j]);
		}
		printf("\n");
	}
	printf("========================================\n");
}

void freeScoreMatrix(mtype **scoreMatrix, int sizeB) {
	int i;
	for (i = 0; i < (sizeB + 1); i++)
		free(scoreMatrix[i]);
	free(scoreMatrix);
}

int main(int argc, char ** argv) {
	// sequence pointers for both sequences
	// seqC armazena todos os caracteres distintos de seqB
	char *seqA, *seqB, *seqC;
	//matriz de suporte para armazernar os valores de K
	//**p_matrix;

	// sizes of both sequences
	int sizeA, sizeB, sizeC;

	//read both sequences
	seqA = read_seq("fileC.in");
	sizeA = strlen(seqA);
	for (int i = 0; i < sizeA ; i++) {
		if ((unsigned char)seqA[i] < 32 || (unsigned char)seqA[i] > 126) {
			printf("OI, Caractere inválido encontrado: seqA[%d] = '%c' (%d)\n", i, seqA[i], (int)seqA[i]);
		}
	}
	seqB = read_seq("fileD.in");
	for (int i = 0; i < sizeB ; i++) {
		if ((unsigned char)seqB[i] < 32 || (unsigned char)seqB[i] > 126) {
			printf("OI, Caractere inválido encontrado: seqB[%d] = '%c' (%d)\n", i, seqB[i], (int)seqB[i]);
		}
	}
	seqC = "ATCG";


	//find out sizes

	sizeB = strlen(seqB);
	sizeC = strlen(seqC);
	printf("%d %d %d\n", sizeA, sizeB, sizeC);

	// allocate LCS score matrix
	mtype ** scoreMatrix = allocateScoreMatrix(sizeB, sizeA);
	mtype ** p_matrix = allocateScoreMatrix(sizeC, sizeB);

	//initialize LCS score matrix
	initScoreMatrix(scoreMatrix, sizeA, sizeB);

	//fill up the rest of the matrix and return final score (element locate at the last line and collumn)
	calcula_p(p_matrix, seqB, sizeB, seqC, sizeC);
	mtype score = LCS_par(scoreMatrix, p_matrix, sizeA, sizeB, sizeC, seqA, seqB, seqC);

	/* if you wish to see the entire score matrix,
	 for debug purposes, define DEBUGMATRIX. */
#ifdef DEBUGMATRIX
	printMatrix(seqA, seqB, scoreMatrix, sizeA, sizeB);
#endif

	//print score
	printf("\nScore: %d\n", score);

	//free score matrix
	freeScoreMatrix(scoreMatrix, sizeB);
	freeScoreMatrix(p_matrix, sizeB);

	return EXIT_SUCCESS;
}
