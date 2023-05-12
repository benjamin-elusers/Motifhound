/* +------------------------------------------------------------+ */
/* |  Date    : 16/05/2011                                      | */
/* |  Author : B. Dubreuil                                      | */
/* |  Contact : dubreuil.benjamin@hotmail.com                   | */
/* |  Centre de Bioinformatique et Genomique Robert Cedergren   | */
/* |  MICHNICK LAB - Université de Montréal (UDEM)              | */
/* +------------------------------------------------------------+ */
/*============================================================================================*/
/*##     ##  #######  ######## #### ######## ##     ##  #######  ##     ## ##    ## ########  */
/*###   ### ##     ##    ##     ##  ##       ##     ## ##     ## ##     ## ###   ## ##     ## */
/*#### #### ##     ##    ##     ##  ##       ##     ## ##     ## ##     ## ####  ## ##     ## */
/*## ### ## ##     ##    ##     ##  ######   ######### ##     ## ##     ## ## ## ## ##     ## */
/*##     ## ##     ##    ##     ##  ##       ##     ## ##     ## ##     ## ##  #### ##     ## */
/*##     ## ##     ##    ##     ##  ##       ##     ## ##     ## ##     ## ##   ### ##     ## */
/*##     ##  #######     ##    #### ##       ##     ##  #######   #######  ##    ## ########  */
/*                               Over Representation Computation                              */
/*============================================================================================*/
/* Takes a Motif Enumeration files and apply the Cumulative HyperGeometric Distribution to obtain for each motif a Pvalue associated.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include "Judy.h"

#define MIN(a,b)   (a < b ? a : b)
#define MAX(a,b)   (a > b ? a : b)
#define MAXLINE 512 /* max string (line) length */
#define MAXHEADER 512 /* max string (line) length */
#define JUDYERROR_SAMPLE 1
#define LINELENGTH 50
#define strdup(x) my_strdup(x)

uint8_t   Index[MAXLINE];		 /* string to insert */

void		my_sprintf(int,char *, long int);												/*	Print a combination (a mask) depending on the size of the motif									*/
void *		safemalloc(size_t); 															/*	malloc Function but check if the memory allocation has not returned NULL pointer				*/
void *		safecalloc(int, size_t); 														/*	calloc Function but check if the memory allocation has not returned NULL pointer				*/
void 		safefree(void *);																/*	free Function but check if the memory deallocation has returned a NULL pointer					*/
char *		my_strdup(const char *);														/*	Save a copy of a string																			*/
int			parse_csv(char *, char *, double, Pvoid_t *, int);					/*	Comma-Separated-Values File Parser																*/
long double	HyperGeometric_Cum_lower(long double, long double, long double, long double);	/*	Compute a p-value from the Hypergeometric cumulative Distribution for the lower tail P(X<=x) )	*/
long double	HyperGeometric_Cum_upper(long double, long double, long double, long double);	/*	Compute a p-value from the Hypergeometric cumulative Distribution for the upper tail P(X>=x) )	*/
int			Compare_Prot_vs_Seqs(double, Pvoid_t, Pvoid_t, char *, int, int); 				/*	Write in a file,for each motif :																
																								|	- the motif pattern,																		|
																								|	- k, his occurrences in the Set File,														|
																								|	- m, his occurrences in the Proteome File,          										| 																									|	- n, the Number of sequences in the Set File,												|
																								|	- N, the Number of sequences in the Proteome File,											|
																								|	- result, the associated pvalue computed from the Hypergeometric cumulative ditribution		|
																								|	- comb, the combination space of a pattern ( equals to 20^(# of non-wildcard positions))	|
																								|	- nb_patt, the number of different pattern with a fixed number of wildcard positions 		|
																								|	- the adjusted pvalue ( pvalue divided or multiplied by 				|	*/


/*******************************************************************************************************/
/**********************************ALLOC/DEALLOCATION OF MEMORY SPACE***********************************/
/*******************************************************************************************************/

void * safecalloc(int taille,size_t s)
/* calloc with test on the memory allocation */
{
	void *mptr = calloc(taille,s);
	if(mptr == NULL) {    exit(ENOMEM);    }
	return mptr;
}

void * safemalloc(size_t s) 
/* malloc with test on the memory allocation*/
{
	void *mptr = malloc(s);
	if(mptr == NULL) {    exit(ENOMEM);    }
	return mptr;
}

void safefree(void *p){
    free(p);
    if(p != NULL){ p=NULL; }
}

char *my_strdup(const char *str) {
/* Copy the Character Array given in argument in a new array with the same size (Character Array Duplication) */ 
	int n = strlen(str) + 1;
	char *dup = safemalloc(n*sizeof(char));
	if(dup) {    strcpy(dup, str);    }
	dup[n-1]='\0';
	return dup;
}

/*******************************************************************************************************/
/*************************CALCULATOR OF CUMULATIVE HYPERGEOMETRIC DISTRIBUTION**************************/
/*******************************************************************************************************/

/*Cumulative hypergeometric distribution*/
long double HyperGeometric_Cum_lower(long double k, long double n, long double m, long double N) {
/* k: successes in the set, n: size of the set, m: successes in the population, N: size of the population */

	long double hg=0.0f;  long double i;  long double j;  long double bc=1.0f; /* long double  bc_tmp = 1.0f; */

	if(m <= N && k <= n && k <= m && n <= (N - m) && (n - k) <= (N - m)) {

		for (i=k; i<= 0; i-=1.0f) {
/*			bc_tmp = bc;*/
			bc = 1.0f;
/*			printf("i=%.0Lf \t",i);*/
			for (j=0.0f; j<n; j+=1.0f) {
				if (j < i)   bc *= (m - i + j + 1.0f) / (j + 1.0f);
				if (j < n-i) bc *= (N-m - n+i + j + 1.0f) / (j + 1.0f);
				if (j < n)   bc /= (N - n + j + 1.0f) / (j + 1.0f);
			}
			hg += bc;
/*			printf("P(X=%.0Lf) = %.80Lf \t P(X>=%.0Lf) = %.80Lf \n",i,bc,i,hg);*/
		}
/*		Return pvalue*/
		return (1-hg);
	}else{
		return 0.0f;
	}
}

/*Cumulative hypergeometric distribution*/
long double HyperGeometric_Cum_upper(long double k, long double n, long double m, long double N) {
/* k: successes in the set, n: size of the set, m: successes in the population, N: size of the population */

	long double hg=0.0f;  long double i;  long double j;  long double bc=1.0f; /* long double  bc_tmp = 1.0f; */

	if(m <= N && k <= n && k <= m && n <= (N - m) && (n - k) <= (N - m)) {

		for (i=k; i<= MIN(n,m); i+=1.0f) {
/*			bc_tmp = bc;*/
			bc = 1.0f;
/*			printf("i=%.0Lf \t",i);*/
			for (j=0.0f; j<n; j+=1.0f) {
				if (j < i)   bc *= (m - i + j + 1.0f) / (j + 1.0f);
				if (j < n-i) bc *= (N-m - n+i + j + 1.0f) / (j + 1.0f);
				if (j < n)   bc /= (N - n + j + 1.0f) / (j + 1.0f);
			}
			hg += bc;
/*			printf("P(X=%.0Lf) = %.80Lf \t P(X>=%.0Lf) = %.80Lf \n",i,bc,i,hg);*/
		}
/*		Return pvalue*/
		return (hg);
	}else{
		return 0.0f;
	}
}

/*Binomial coefficient for large numbers*/
long double coeff_binom(long double n,long double k){
	long double r = 1.0f;
	long double d = n - k;
	/* choose the smaller of k and n - k */
	if (d > k) { k = d; d = n - k; }

	while (n > k) {
		r *= n--;
		/* divide (n - k)! as soon as we can to delay overflows */
		while (d > 1.0f && !fmod(r,d)) r /= d--;
	}
	return r;
}


/*******************************************************************************************************/
/*******************************************************************************************************/
/*******************************************************************************************************/

int Count_Wildcards(char * mot, double length){
	
	int p=0; int nb_wild=0;
	
	for(p=1;p<length-1;p++){
		if(mot[p] == '.'){ nb_wild++; }
	}
	
	return (nb_wild);
}


int Compare_Prot_vs_Seqs(double S,Pvoid_t k_Array, Pvoid_t m_Array,char *path_output,int nSeq,int NSeq) {

	int w=0; 
	FILE* output = (FILE*)NULL;
/*		FILE* test = (FILE*)NULL;*/

	/*int inutile =0;*/
	long double result = 0.0f;
/*	long double result_adj = 0.0f;*/
/*	long double comb = 0.0f;*/
/*	long double nb_patt = 0.0f;*/
	PWord_t   k,m;     /* Judy array element. */

	output = fopen(path_output, "a");
	if (!output)	{	fprintf(stderr,"Error! Unable to open the results file \"%s\".\n",path_output);	exit(1);	}
	
/*	test = fopen("test.txt", "a");*/
/*	if (!test)	{	fprintf(stderr,"Error! Unable to open the test file \"test.txt\".\n");	exit(1);	}*/

	/* Prints the motifs sorted lexicographically */
	Index[0] = '\0';
	JSLF(k, k_Array, Index);	/* get first string */
/*	printf("%s\n",(char *)(Index));*/
	while (k != NULL) {
		/* w=Count_Wildcards( (char*)Index, S);*/
		/* comb=pow(20,S-w);  nb_patt=coeff_binom(S-2,w);*/
		JSLG(m, m_Array, Index);
		if(m != NULL) {
			if((int)(*m) >= (int)(*k)){
				result = HyperGeometric_Cum_upper((long double)(*k),(long double)(nSeq),(long double)(*m),(long double)(NSeq));
				fprintf(output,"%s\t%d\t%d\t%d\t%d\t%.8Le\n", Index,(int)(*k),(int)(*m),(int)(nSeq),NSeq,result);
/*				result_adj = result * nb_patt;*/
/*				fprintf(output,"%s\t%d\t%d\t%d\t%d\t%.8Le\t%.2Le\t%.8Le\t%d\n", Index,(int)(*k),(int)(*m),(int)(nSeq),NSeq,result,nb_patt,result_adj,w);*/
/*				fprintf(test,"m>=k:%s,%d,%d,%d,%d,%.8Le\n",(char *)(Index),(int)(*k),(int)(*m),nSeq,NSeq,result);*/
			}else{
				result = HyperGeometric_Cum_upper((long double)(*k),(long double)(nSeq),(long double)(*m),(long double)(NSeq));
/*				result_adj = result * nb_patt;*/
				fprintf(output,"%s\t%d\t%d\t%d\t%d\t%.8Le | ( WARNING: This motif occurs more in the Set Sequences than in the Background Sequences )\n", Index,(int)(*k),(int)(*m),(int)(nSeq),NSeq,result);
/*				fprintf(output,"%s\t%d\t%d\t%d\t%d\t%.8Le\t%.2Le\t%.8Le\t%d | ( WARNING: This motif occurs more in the Set Sequences than in the Background Sequences )\n", Index,(int)(*k),(int)(*m),(int)(nSeq),NSeq,result,nb_patt,result_adj,w);*/
/*				fprintf(test,"m<k: %s,%d,%d,%d,%d,%.8Le\n",(char *)(Index),(int)(*k),(int)(*m),nSeq,NSeq,result);*/
			}
		}else{
			result = HyperGeometric_Cum_upper((long double)(*k),(long double)(nSeq),(long double)(0+*k),(long double)(NSeq));
/*			result_adj = result * nb_patt;*/
			fprintf(output,"%s\t%d\t%d\t%d\t%d\t%.8Le | ( WARNING: This motif might occur less than 3 times in the Background Sequences )\n", Index,(int)(*k),(int)(*m),(int)(nSeq),NSeq,result);
/*			fprintf(output,"%s\t%d\t%d\t%d\t%d\t%.8Le\t%.2Le\t%.8Le\t%d | ( WARNING: This motif might occur less than 3 times in the Background Sequences )\n", Index,(int)(*k),(int)(0+*k),(int)(nSeq),NSeq,result,nb_patt,result_adj,w);*/
/*			fprintf(test,"m=0: %s,%d,%d,%d,%d,%.8Le\n",(char *)(Index),(int)(*k),(int)(*m),nSeq,NSeq,result);*/
		}
		JSLN(k, k_Array, Index);   /* get next string */
	}
	fclose(output);
/*	fclose(test);*/
	return 0;
}

int parse_csv(char *path_csvfile, char *separateur,double s, Pvoid_t * success_Array, int OnevsAll) {

	FILE* csvfile = (FILE*)NULL;
	char line [MAXLINE]; char *pcourant=NULL; /*int irand=0;*/
	PWord_t success;     /* Judy array element. */
/*	int inutile = 0;*/
/*		FILE* test = (FILE*)NULL;*/
/*		char testfile[256];*/
	csvfile = fopen(path_csvfile, "r");
	/* Test si le fichier peut etre ouvert en lecture */
	if (!csvfile)	{	fprintf(stderr,"Error! Unable to open the CSV formatted file of Enumeration \"%s\".\n",path_csvfile);	exit(1);	}

	/* initialize random seed: */
/*	srand ( time(NULL) );*/
	/* generate random number: */
/*	irand = rand() % 10 + 1;*/
/*	sprintf(testfile,"csv%d.txt",irand);*/
/*	test = fopen(testfile, "a");*/
/*	if (!test)	{	fprintf(stderr,"Error! Unable to open the test file \"%s\".\n",testfile);	exit(1);	}*/

/*	printf("%s",path_csvfile);*/

	while(fgets(line,MAXHEADER,csvfile)) {
		if( line[0] != '#' ){
	/*		printf("%s",line);*/

			pcourant = strtok(line, separateur);
	/*		printf("pcourant=%s \t",pcourant);*/
			if (pcourant != NULL){
				strcpy((char *)(Index),pcourant); Index[(int)(s)]='\0';
			}
	/*		printf("Motif=%s \t",(char *)(Index));*/

			/********************************/
			/* All occurrences per sequence */
			/********************************/
			pcourant = strtok(NULL, separateur);
			if(OnevsAll == 0){
				if (pcourant != NULL){
					JSLI(success, *(success_Array),Index);   /* store string into array */
					(*success)=atoi(pcourant);
				}
			}
	/*		printf("pcourant=%s \t",pcourant);*/
	/*		printf("pcourant=%s \t",pcourant);*/
	
			/*******************************/
			/* One occurrence per sequence */
			/*******************************/
			pcourant = strtok(NULL, "\n");
			if(OnevsAll != 0){
				if (pcourant != NULL){
					JSLI(success, *(success_Array),Index);   /* store string into array */
					(*success)=atoi(pcourant);
				}
			}
/*			fprintf(test,"%s,%d\n",Index,(int)(*success));*/
	/*		printf("\n");*/

		}
	}
/*	scanf("%d",&inutile);*/
	fclose(csvfile);
	return 0;
}

int main (int argc,char **argv){

	/* DECLARATION */
	int nSeq; int NSeq; int OnevsAll; /* int i; int ii;	int inutile; */

	double size; 

	char * Prot_csvfile = NULL; char * Seqs_csvfile = NULL; char * Output_path = NULL;

	Pvoid_t   P_mArray = (PWord_t)NULL;  /* Judy array. */
	Pvoid_t   S_kArray = (PWord_t)NULL;  /* Judy array. */

/*	PWord_t   k,m;     Judy array element. */

	Word_t    Rc_word;		     /* size of JudySL array. */

	/* INITIALIZATION */
	size=0.0f;
	nSeq=0; NSeq=0; OnevsAll=0;/* i=0; ii=0; inutile = 0; */

	/* AFFECTATION */
	Prot_csvfile =  strdup(argv[1]);
	Seqs_csvfile =  strdup(argv[2]); 
	Output_path  =  strdup(argv[3]); 
	size = atof(argv[4]);
	nSeq = atoi(argv[5]);
	NSeq = atoi(argv[6]);
	if(argc == 8){
		OnevsAll = atoi(argv[7]);
	}


/*	printf("\n");*/

/*	printf("size %.0f\t",size);*/
/*	printf("it %d\n",it);*/
/*	printf("INPUT :\n Sequences => %s \n Proteome => %s \n",Seqs_csvfile,Prot_csvfile);*/

/*	scanf("%d",&inutile);*/
/*	printf(" #Seq %d\t\n",nSeq);*/

	parse_csv(Prot_csvfile,"\t",size,&P_mArray,OnevsAll);
	parse_csv(Seqs_csvfile,"\t",size,&S_kArray,OnevsAll);

	/* Calculate the hypergeometric cumulative distribution for each motif -> Pvalue */
	Compare_Prot_vs_Seqs(size,S_kArray,P_mArray,Output_path,nSeq,NSeq);
	JSLFA(Rc_word,S_kArray);
	JSLFA(Rc_word,P_mArray);


	safefree(Prot_csvfile);
	safefree(Seqs_csvfile); 
	safefree(Output_path); 

/*	printf("\n");*/
	return EXIT_SUCCESS;
}

