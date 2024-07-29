/*
*********************************************************************************************************
*	About: The code implements Quantum First - Last Filtering based Multiple Pattern Matching (QFL_MPM)	*
*	Usage: Run in command prompt "QFL_MPM.exe <input_file>"												*
*	1. 	<input_file> contains character subset of biological gene text sequence of SARS-CoV-2 for Human	*
* 		File sizes {128, 256, 512} characters are intentioally prepared for feasible QuEST Simulation 	*
*	2.	QFL filtering and MPM searching is performed for multiple patterns P1 and P2 of length {5, 4}	*
		Multiple Quantum Core (QCore = 2) realization based QFL filtering and MPM searching simulation	*
*		Quantum circuits are designed for search patterns P1 = "A C T A G" & P2 = "G T T A" within text	*
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

//Binary Adder
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
	
	//STARTS - Reading and Storing Input File in the array "T"
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
	int T_size = count;
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
	
	//The arbitrary patterns P1 = "ACTAG" & P2 = "GTTA" used for searching within the text "T"
	const int P1_size = 5, P2_size = 4;
	char P1[] =  { 'A', 'C', 'T', 'A', 'G' };
	char P2[] =  { 'G', 'T', 'T', 'A'};

	printf("\nDNA sequence is:\n");
	for(int i = 1; i <= count; i++)
	{
		printf("%c", T[i - 1]);
		if(i % 32 == 0) printf("\n");
	}
	printf("\nPatterns to search are: ");
	for(int i = 0; i < P1_size; i++)
		printf("%c", P1[i]);
	printf(", ");
	for(int i = 0; i < P2_size; i++)
		printf("%c", P2[i]);
	printf("\n");

	//STARTS - Preliminaries for quantum circuit corresponding to unitaries U_QFil, U_TPos and U_QAdd
	//NOTE: This is not the efficient approach, however, can be improved for unrestricted simulation
	
	int n = (int)(log(count)/log(2));

	int **anf1_LA, *anf_size1_LA;
	int **anf2_LA, *anf_size2_LA;
	int *f1_LA = (int *)malloc(sizeof(int) * (int)pow(2, n));
	int *f2_LA = (int *)malloc(sizeof(int) * (int)pow(2, n));
    for(int i = 0; i < count; i++) //count = 2^n
    {
		if(tolower(T[i]) == tolower(P1[0]) && tolower(T[i + P1_size - 1]) == tolower(P1[P1_size - 1]))
			f1_LA[i] = bin_add(i + P1_size - 1, _2s_complement(P1_size - 1, n), n);
		else
			f1_LA[i] = bin_add(i, _2s_complement(i, n), n);
		
		if(tolower(T[i]) == tolower(P2[0]) && tolower(T[i + P2_size - 1]) == tolower(P2[P2_size - 1]))
			f2_LA[i] = bin_add(i + P2_size - 1, _2s_complement(P2_size - 1, n), n);
		else
			f2_LA[i] = bin_add(i, _2s_complement(i, n), n);
    }
	anf_calc(f1_LA, n, &anf1_LA, &anf_size1_LA);
	anf_calc(f2_LA, n, &anf2_LA, &anf_size2_LA);
	free(f1_LA);
	free(f2_LA);

	//ENDS - Preliminaries for quantum circuit corresponding to unitaries U_QFil, U_TPos and U_QAdd


	//STARTS - Quantum First - Last Filtering (QFL) Call based on Multiple Quantum Core (QCore = 2)
	//Quantum Environment (env) Realizing Two Quantum System (qubits1 and qubits2) = Simulates Two Quantum Core over Multithreaded Classical Core

	int res_count1 = 0, res_count2 = 0;
	int realloc_size1 = 1, realloc_size2 = 1;
	int *Loc1 = (int *)malloc(sizeof(int) * realloc_size1);
	int *Loc2 = (int *)malloc(sizeof(int) * realloc_size2);

	QuESTEnv env = createQuESTEnv();
	Qureg qubits1 = createQureg(2 * n, env);
	Qureg qubits2 = createQureg(2 * n, env);
	printf("\nQuantum Parameters for Filtering Part (same for all cores):\n");
	reportQuregParams(qubits1);
	reportQuESTEnv(env);
	printf("\n");

	//Single qubit unitary for XOR gate
	ComplexMatrix2 ux = {
		.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
	};

	for(int repeat = 0; repeat < 500; repeat++)
	{
		initZeroState(qubits1);
		initZeroState(qubits2);

		for(int i = n; i < 2 * n; i++)
		{
			hadamard(qubits1, i);
			hadamard(qubits2, i);
		}

		//Evaluation of U_TPos and U_QAdd by below quantum circuit over Simulated QCore 1
		for(idx = 0; idx < n; idx++)
		{
			for(int i = 0; i < anf_size1_LA[idx]; i++)
			{
				if(anf1_LA[idx][i] == 0)
				{
					pauliX(qubits1, idx);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * n);
				int ctrl_size = 0;
				int term = anf1_LA[idx][i];
				int qb = 0;
				while(term)
				{
					if(term & 1) ctrls[ctrl_size++] = (qb + n);
					term >>= 1;
					qb++;
				}
				multiControlledUnitary(qubits1, ctrls, ctrl_size, idx, ux);
				free(ctrls);
			}
		}

		//Evaluation of U_TPos and U_QAdd by below quantum circuit over Simulated QCore 2
		for(idx = 0; idx < n; idx++)
		{
			for(int i = 0; i < anf_size2_LA[idx]; i++)
			{
				if(anf2_LA[idx][i] == 0)
				{
					pauliX(qubits2, idx);
					continue;
				}
				int *ctrls = (int *)malloc(sizeof(int) * n);
				int ctrl_size = 0;
				int term = anf2_LA[idx][i];
				int qb = 0;
				while(term)
				{
					if(term & 1) ctrls[ctrl_size++] = (qb + n);
					term >>= 1;
					qb++;
				}
				multiControlledUnitary(qubits2, ctrls, ctrl_size, idx, ux);
				free(ctrls);
			}
		}

		for(int i = n; i < 2 * n; i++)
		{
			hadamard(qubits1, i);
			hadamard(qubits2, i);
		}

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
				prob = getProbAmp(qubits1, i);
				if(prob != 0.0)
					printf("\ni = %d %d has prob %f", (i >> n), (i & mask), prob);
				prob = getProbAmp(qubits2, i);
				if(prob != 0.0)
					printf(" | %d %d has prob %f", (i >> n), (i & mask), prob);
			}
		}
		*/

		//Counting and storing likely indexes
		/*
		qreal prob;
		for(int i = 0; i < (int)pow(2, n); i++)
		{
			prob = getProbAmp(qubits1, i);
			if(prob != 0.0)
			{
				if(realloc_size1 == res_count1)
				{
					realloc_size1 *= 2;
					Loc1 = (int *)realloc(Loc1, sizeof(int) * realloc_size1);
				}
				Loc1[res_count1++] = i;
			}
			prob = getProbAmp(qubits2, i);
			if(prob != 0.0)
			{
				if(realloc_size2 == res_count2)
				{
					realloc_size2 *= 2;
					Loc2 = (int *)realloc(Loc2, sizeof(int) * realloc_size2);
				}
				Loc2[res_count2++] = i;
			}
		}
		*/

		///*
		int idx1 = 0, idx2 = 0;
		for(int i = 0; i < n; i++)
		{
			int outcome = measure(qubits1, i);
			if(outcome) idx1 ^= (outcome << i);
			outcome = measure(qubits2, i);
			if(outcome) idx2 ^= (outcome << i);
		}
		printf("#%d: Measuring last 'n' qubits return - %d | %d\n", repeat + 1, idx1, idx2);
		if(!in(Loc1, res_count1, idx1))
		{
			if(realloc_size1 == res_count1)
			{
				realloc_size1 *= 2;
				Loc1 = (int *)realloc(Loc1, sizeof(int) * realloc_size1);
			}
			Loc1[res_count1++] = idx1;
		}
		if(!in(Loc2, res_count2, idx2))
		{
			if(realloc_size2 == res_count2)
			{
				realloc_size2 *= 2;
				Loc2 = (int *)realloc(Loc2, sizeof(int) * realloc_size2);
			}
			Loc2[res_count2++] = idx2;
		}
		//*/
	}


	//Free Memory by Algebraic Normal Form (ANF) for Both Simulated QCore 1 and QCore2
	for(int i = 0; i < n; i++)
		free(anf1_LA[i]);
	free(anf1_LA);
	for(int i = 0; i < n; i++)
		free(anf2_LA[i]);
	free(anf2_LA);

	destroyQureg(qubits1, env);
	destroyQureg(qubits2, env);
	destroyQuESTEnv(env);

	//ENDS - Quantum First - Last Filtering (QFL) Method Call
	//Showing likely filtered (non-trivial along with a trivial) indices of Simulated QCore 1

	printf("\nLikely indexes are (from core 1):\n");
	for(int i = 0; i < res_count1; i++)
		printf("Index %d\n", Loc1[i]);

	//Showing likely filtered (non-trivial along with a trivial) indices of Simulated QCore 2	
	printf("\nLikely indexes are (from core 2):\n");
	for(int i = 0; i < res_count2; i++)
		printf("Index %d\n", Loc2[i]);



	//STARTS - Multiple Pattern Matching (MPM) based on Multiple Quantum Core (QCore = 2) through Exact - Parallel Match U_Comp & Grover's Search Operator (GSO) Call
	//Quantum Environment (env) Realizing Two Quantum System (qubits1 and qubits2) = Simulates Two Quantum Core over Multithreaded Classical Core

    if(res_count1 < 2)
    {
    	printf("Grover's Searching is not required since only %d index is found!\n", res_count1);
    	if(res_count1 == 1 && is_at(T, count, Loc1[0], P1, P1_size)) printf("Classically checked that the pattern is at index %d.\n", Loc1[0]);
    	else printf("The pattern is not at index %d.\n", Loc1[0]);
    }
	else
	{
		//STARTS - Quantum Grover's Search Operator (GSO) Call for Simulated QCore 1

		env = createQuESTEnv();

		int req_qbs1 = ceil(log(res_count1)/log(2));
		qubits1 = createQureg(req_qbs1, env);
		initZeroState(qubits1);
		printf("\nQuantum Parameters for Grover's Part on core 1:\n");
		reportQuregParams(qubits1);
		reportQuESTEnv(env);

		count = 0;
		ComplexMatrixN e1 = createComplexMatrixN(req_qbs1);
		for(int i = 0; i < (int)pow(2, req_qbs1); i++)
		{
			if(i < res_count1 && is_at(T, T_size, Loc1[i], P1, P1_size))
			{
				e1.real[i][i] = -1;
				count++;
			}
			else e1.real[i][i] = 1;
		}

		int *targs1 = (int *)malloc(sizeof(int) * n);
		for(int i = 0; i < req_qbs1; i++)
			targs1[i] = i;
		int times1 =  ceil(3.14 * (pow(2, req_qbs1 / 2) / sqrt(count)) / 4);
		printf("\nRunning Grover's Algorithm on filtered indexes (core 1)..\n");
		for(int i = 0; i < req_qbs1; i++)
			hadamard(qubits1, i);
		for(int gi = 0; gi < times1; gi++)
		{
			//Marking
			multiQubitUnitary(qubits1, targs1, req_qbs1, e1);
			
			//Diffusion
			for(int i = 0; i < req_qbs1; i++)
				hadamard(qubits1, i);
			for(int i = 0; i < req_qbs1; i++)
				pauliX(qubits1, i);
			multiControlledPhaseFlip(qubits1, targs1, req_qbs1);
			for(int i = 0; i < req_qbs1; i++)
				pauliX(qubits1, i);
			for(int i = 0; i < req_qbs1; i++)
				hadamard(qubits1, i);
		}

		qreal max = 0.0;
		qreal prob = 0.0;
		for(int i = 0; i < (int)pow(2, req_qbs1); i++)
		{
			prob = getProbAmp(qubits1, i);
			if(max <= prob) max = prob;
		}
		for(int i = 0; i < (int)pow(2, req_qbs1); i++)
		{
			prob = getProbAmp(qubits1, i);
			if(fabs(max - prob) <= 0.0000001) printf("Correct Index %d\n", Loc1[i]);
		}

		destroyComplexMatrixN(e1);
		destroyQureg(qubits1, env);
		destroyQuESTEnv(env);

		//ENDS - Quantum Grover's call
	}

	if(res_count2 < 2)
    {
    	printf("Grover's Searching is not required since only %d index is found!\n", res_count2);
    	if(res_count2 == 1 && is_at(T, count, Loc2[0], P2, P2_size)) printf("Classically checked that the pattern is at index %d.\n", Loc2[0]);
    	else printf("The pattern is not at index %d.\n", Loc2[0]);
    }
	else
	{
		//STARTS - Quantum Grover's Search Operator (GSO) Call for Simulated QCore 2

		env = createQuESTEnv();

		int req_qbs2 = ceil(log(res_count2)/log(2));
		qubits2 = createQureg(req_qbs2, env);
		initZeroState(qubits2);
		printf("\nQuantum Parameters for Grover's Part on core 2:\n");
		reportQuregParams(qubits2);
		reportQuESTEnv(env);

		count = 0;
		ComplexMatrixN e2 = createComplexMatrixN(req_qbs2);
		for(int i = 0; i < (int)pow(2, req_qbs2); i++)
		{
			if(i < res_count2 && is_at(T, T_size, Loc2[i], P2, P2_size))
			{
				e2.real[i][i] = -1;
				count++;
			}
			else e2.real[i][i] = 1;
		}

		int *targs2 = (int *)malloc(sizeof(int) * n);
		for(int i = 0; i < req_qbs2; i++)
			targs2[i] = i;
		int times2 =  ceil(3.14 * (pow(2, req_qbs2 / 2) / sqrt(count)) / 4);
		printf("\nRunning Grover's Algorithm on filtered indexes (core 2)..\n");
		for(int i = 0; i < req_qbs2; i++)
			hadamard(qubits2, i);
		for(int gi = 0; gi < times2; gi++)
		{
			//Marking
			multiQubitUnitary(qubits2, targs2, req_qbs2, e2);
			
			//Diffusion
			for(int i = 0; i < req_qbs2; i++)
				hadamard(qubits2, i);
			for(int i = 0; i < req_qbs2; i++)
				pauliX(qubits2, i);
			multiControlledPhaseFlip(qubits2, targs2, req_qbs2);
			for(int i = 0; i < req_qbs2; i++)
				pauliX(qubits2, i);
			for(int i = 0; i < req_qbs2; i++)
				hadamard(qubits2, i);
		}

		qreal max = 0.0;
		qreal prob = 0.0;
		for(int i = 0; i < (int)pow(2, req_qbs2); i++)
		{
			prob = getProbAmp(qubits2, i);
			if(max <= prob) max = prob;
		}
		for(int i = 0; i < (int)pow(2, req_qbs2); i++)
		{
			prob = getProbAmp(qubits2, i);
			if(fabs(max - prob) <= 0.0000001) printf("Correct Index %d\n", Loc2[i]);
		}

		destroyComplexMatrixN(e2);
		destroyQureg(qubits2, env);
		destroyQuESTEnv(env);

	//ENDS - Multiple Pattern Matching (MPM) based on Multiple Quantum Core (QCore = 2) through Exact - Parallel Match U_Comp & Grover's Search Operator (GSO) Call
	}

    return 0;
}