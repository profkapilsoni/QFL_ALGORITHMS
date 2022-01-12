/*
*********************************************************************************************************
*	About: The code implements Quantum First - Last Filtering based Single Pattern Matching (QFL_SPM)	*
*	Usage: Run in command prompt "QFL_SPM.exe <input_file>"												*
*	1. 	<input_file> contains character subset of biological gene text sequence of SARS-CoV-2 for Human	*
* 		File sizes {128, 256, 512} characters are intentioally prepared for feasible QuEST Simulation 	*
*	2.	QFL based filtering and SPM searching is performed for the fixed pattern P of the length 5		*
*		Quantum algorithmic circuit is designed for a searching pattern P = "A C T A G" within text		*
*********************************************************************************************************
*/

// Instead of applying unitary circuit we use their corresponding unitary matrices
// Some realization of quantum circuits are done classically, instead of quantum
// NOTE: This is done to make the QuEST simulation as memory efficient and effective

#include <stdio.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

//Check if an integer element is in the given array: Returns 1 if the element is present and 0, otherwise
int in(int *a, int size, int itm)
{
	for(int i = 0; i < size; i++)
		if(a[i] == itm) return 1;
	return 0;
} 


//Used in building Algebraic Normal Form (ANF) 
int check_term(int t, int a, int width)
{
	for (int i = 0; i < width; i++)
	{
		if ((t & (1 << i)) == 0 && (a & (1 << i)) != 0) return 0;
	}
	return 1;
}

//Build / Construct Algebraic Normal Form (ANF): Realized for Quantum Memory, Adder Operation, and Exact - Parallel Match
void anf_calc(int f[], int n, int ***anf, int **anf_size)
{
	*anf = (int **)malloc(sizeof(int *) * n);
	for(int i = 0; i < n; i++)
		(*anf)[i] = (int *)malloc(sizeof(int) * (int)pow(2, n));
	*anf_size = (int *)malloc(sizeof(int) * n);
	for(int i = 0; i < n; i++)
		(*anf_size)[i] = 0;

	for (int i = n - 1; i >= 0; i--)
	{
		for (int term = 0; term < (int)pow(2, n); term++)
		{
			int count = 0;
			for (int fi = 0; fi < (int)pow(2, n); fi++)
			{
				if (check_term(term, fi, n) && (f[fi] & (1 << i)) != 0) count++;
			}
			if (count % 2 == 1) (*anf)[i][(*anf_size)[i]++] = term;
		}
	}
}

//Return 2's complement of an integer of "n" bits
int _2s_complement(int num, int bits)
{
	for(int i = 0; i < bits; i++)
		num ^= (1 << i);
	return (num + 1);
}

//Binary Adder Operation
int bin_add(int a, int b, int width)
{
	int carry = 0, sum = 0;
	for(int i = 0; i < width; i++)
	{
		int _a = ((a & (1 << i)) >> i);
		int _b = ((b & (1 << i)) >> i);
		sum ^= ((_a ^ _b ^ carry) << i);
		if((_a + _b + carry) >= 2)
			carry = 1;
		else carry = 0;
	}
	return sum;
}

int is_at(char *T, int T_size, int idx, char *P, int P_size)
{
	int found = 1;
	int temp_i = 0;
	if((idx + P_size) > T_size) return 0;
	for(int i = idx; i < idx + P_size; i++)
		if(toupper(T[i]) != P[temp_i++]) found = 0;
	return found;
}

int main (int narg, char *varg[])
{
	if(narg != 2)
	{
		printf("Usage: %s <input_file>", varg[0]);
		return 0;
	}
	
	//STARTS - Reading and Storing the Input File in the text array "T"
	//The code assumes that the number of elements in the input file is in power of 2
	//for the purpose of easier processing and binary equivalent matching

	FILE *fp = fopen(varg[1], "r");
	int count = 0;
	char ch;
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			count++;
	}
	fclose(fp);
	
	char *T = (char *)malloc(sizeof(char) * count);
	int idx = 0;
	fp = fopen(varg[1], "r");
	while((ch = fgetc(fp)) != EOF)
	{
		if(tolower(ch) == 'a' || tolower(ch) == 't' || tolower(ch) == 'g' || tolower(ch) == 'c')
			T[idx++] = ch;
	}
	fclose(fp);

	//ENDS - Reading and Storing Input File in the text array "T"
	
	//An arbitrary pattern P = "ACTAG" used for searching within the text "T"
	const int P_size = 5;
	char P[] =  { 'A', 'C', 'T', 'A', 'G' };

	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= count; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}
	printf("Pattern to search: ");
	for(int i = 0; i < P_size; i++)
		printf("%c", P[i]);
	printf("\n");

	//STARTS - Preliminaries for quantum circuit corresponding to unitaries U_QFil, U_TPos and U_QAdd
	//NOTE: This is not the efficient approach, however, can be improved for unrestricted simulation
	
	int n = (int)(log(count)/log(2));

	//int **anf_FI, *anf_size_FI;
	int **anf_LA, *anf_size_LA;
	int *f_LA = (int *)malloc(sizeof(int) * (int)pow(2, n));
    for(int i = 0; i < count; i++) //count = 2^n
    {
		if(tolower(T[i]) == tolower(P[0]) && tolower(T[i + P_size - 1]) == tolower(P[P_size - 1]))
			f_LA[i] = bin_add(i + P_size - 1, _2s_complement(P_size - 1, n), n);
		else
			f_LA[i] = bin_add(i, _2s_complement(i, n), n);
    }
	anf_calc(f_LA, n, &anf_LA, &anf_size_LA);
	free(f_LA);

	//ENDS - Preliminaries for quantum circuit corresponding to unitaries U_QFil, U_TPos and U_QAdd


	//STARTS - Quantum First - Last Filtering (QFL) Call
	//Quantum Environment (env) Realizing Single Quantum System (qubits) = Simulates One Quantum Core over Multithreaded Classical Core

	int res_count = 0;
	int realloc_size = 1;
	int *Loc = (int *)malloc(sizeof(int) * realloc_size);

	QuESTEnv env = createQuESTEnv();
	Qureg qubits = createQureg(2 * n, env);
	printf("\nQuantum Parameters for Filtering Part:\n");
	reportQuregParams(qubits);
	reportQuESTEnv(env);

	int repeat;
	for(repeat = 0; repeat < 500; repeat++)
	{
		initZeroState(qubits);

		//Single qubit unitary for XOR gate
		ComplexMatrix2 ux = {
			.real={{0,1},{1,0}},
			.imag={{0,0},{0,0}}
		};

		for(int i = n; i < 2 * n; i++)
			hadamard(qubits, i);

		//Evaluation of U_TPos and U_QAdd by below quantum circuit
		for(idx = 0; idx < n; idx++)
		{
			for(int i = 0; i < anf_size_LA[idx]; i++)
			{
				if(anf_LA[idx][i] == 0)
				{
					pauliX(qubits, idx);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * n);
				int ctrl_size = 0;
				int term = anf_LA[idx][i];
				int qb = 0;
				while(term)
				{
					if(term & 1) ctrls[ctrl_size++] = (qb + n);
					term >>= 1;
					qb++;
				}
				multiControlledUnitary(qubits, ctrls, ctrl_size, idx, ux);
				free(ctrls);
			}
		}

		for(int i = n; i < 2 * n; i++)
			hadamard(qubits, i);

		//Prints the quantum state before applying Grover's algorithm
		/*
		printf("Required Quantum State is constructed.\n");
		printf("Do you want to view the constructed quantum state?(y/n):");
		ch = getch();
		if(ch == 'y')
		{
			printf("\nConstructed quantum state before applying Grover's search is:");
			qreal prob;
			int mask = (int)pow(2, n) - 1;
			for(int i = 0; i < (int)pow(2, 2 * n); i++)
			{
				prob = getProbAmp(qubits, i);
				if(prob != 0.0)
					printf("\ni = %d %d has prob %f", (i >> n), (i & mask), prob);
			}
		}
		*/

		//Counting and storing likely indexes
		/*
		qreal prob;
		for(int i = 0; i < (int)pow(2, n); i++)
		{
			prob = getProbAmp(qubits, i);
			if(prob != 0.0)
			{
				if(realloc_size == res_count)
				{
					realloc_size *= 2;
					Loc = (int *)realloc(Loc, sizeof(int) * realloc_size);
				}
				Loc[res_count++] = i;
			}
		}
		*/

		///*
		idx = 0;
		for(int i = 0; i < n; i++)
		{
			int outcome = measure(qubits, i);
			if(outcome) idx ^= (outcome << i);
		}
		printf("#%d: Measuring last 'n' qubits return - %d\n", repeat + 1, idx);
		if(!in(Loc, res_count, idx))
		{
			if(realloc_size == res_count)
			{
				realloc_size *= 2;
				Loc = (int *)realloc(Loc, sizeof(int) * realloc_size);
			}
			Loc[res_count++] = idx;
		}
		//*/
	}

	//Free Memory by Algebraic Normal Form (ANF)
	for(int i = 0; i < n; i++)
		free(anf_LA[i]);
	free(anf_LA);

	//ENDS - Quantum First - Last Filtering (QFL) Method Call
	//Showing likely filtered (non-trivial along with a trivial) indices
	
	printf("\nLikely indexes are:\n");
	for(int i = 0; i < res_count; i++)
		printf("Index %d\n", Loc[i]);

	//STARTS - Single Pattern Matching (SPM) through Exact - Parallel Match U_Comp & Grover's Search Operator (GSO) Call
	//Quantum Environment (env) Realizing Single Quantum System (qubits) = Simulates One Quantum Core over Multithreaded Classical Core

    if(res_count < 2)
    {
    	printf("Grover's Searching is not required since only %d index is found!\n", res_count);
    	if(res_count == 1 && is_at(T, count, Loc[0], P, P_size)) printf("Classical checked that the pattern is at index %d.\n", Loc[0]);
    	else printf("The pattern is not at index %d.\n", Loc[0]);
    	return 0;
    }  

    env = createQuESTEnv();

    int req_qbs = ceil(log(res_count)/log(2));
    qubits = createQureg(req_qbs, env);
    initZeroState(qubits);
    printf("\nQuantum Parameters for Grover's Part:\n");
    reportQuregParams(qubits);
    reportQuESTEnv(env);

	int T_size = count;
	count = 0;
	ComplexMatrixN e = createComplexMatrixN(req_qbs);
    for(int i = 0; i < (int)pow(2, req_qbs); i++)
    {
    	if(i < res_count && is_at(T, T_size, Loc[i], P, P_size))
    	{
    		e.real[i][i] = -1;
    		count++;
    	}
    	else e.real[i][i] = 1;
    }

	int *targs = (int *)malloc(sizeof(int) * n);
	for(int i = 0; i < req_qbs; i++)
		targs[i] = i;
    int times =  ceil(3.14 * (pow(2, req_qbs / 2) / sqrt(count)) / 4);

    printf("\nRunning Grover's Algorithm on filtered indexes..\n");
    for(int i = 0; i < req_qbs; i++)
        hadamard(qubits, i);
    for(int gi = 0; gi < times; gi++)
    {
    	//Marking
    	multiQubitUnitary(qubits, targs, req_qbs, e);
    	
        //Diffusion
        for(int i = 0; i < req_qbs; i++)
            hadamard(qubits, i);
        for(int i = 0; i < req_qbs; i++)
            pauliX(qubits, i);
        multiControlledPhaseFlip(qubits, targs, req_qbs);
        for(int i = 0; i < req_qbs; i++)
            pauliX(qubits, i);
        for(int i = 0; i < req_qbs; i++)
            hadamard(qubits, i);
    }

    qreal max = 0.0;
    qreal prob = 0.0;
    for(int i = 0; i < (int)pow(2, req_qbs); i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(max <= prob) max = prob;
    }
    for(int i = 0; i < (int)pow(2, req_qbs); i++)
    {
    	prob = getProbAmp(qubits, i);
    	if(fabs(max - prob) <= 0.0000001) printf("Correct Index %d\n", Loc[i]);
    }

    destroyQureg(qubits, env);
	destroyComplexMatrixN(e);
    destroyQuESTEnv(env);

	//ENDS - Single Pattern Matching (SPM) through Exact - Parallel Match U_Comp & Grover's Search Operator (GSO) Call

    return 0;
}
// END OF CODE
