/* +------------------------------------------------------------+ */
/* |  Date    : 10/05/2012                                      | */
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
/*                                    Motif Scan in Proteome                                  */
/*============================================================================================*/
/* Takes a Motif Enumeration file in input and scan for motifs in a set of Sequences*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include "Judy.h"
#define JUDYERROR_SAMPLE 1
#define MIN2(a,b)   (a < b ? a : b)
#define MAXLINE 512 /* max string (line) length */
#define MAXSEQ 40000 /* max sequence length */
#define MAXHEADER 512 /* max header length */
#define strdup(x) my_strdup(x)

typedef struct fastafile_s {
  FILE *fp;
  char  buffer[MAXSEQ];
} FASTAFILE;

typedef struct {
	char *header;
	char *sequence;
	int Length;
} Fasta_seq;
typedef Fasta_seq * FastaSequences;


void		my_sprintf(int,char *,long int);												/*	Print a combination (a mask) depending on the size of the motif 											*/
void *		safemalloc(size_t); 															/*	malloc Function but check if the memory allocation has not returned NULL pointer							*/
void *		safecalloc(int, size_t); 														/*	calloc Function but check if the memory allocation has not returned NULL pointer							*/
void *		saferealloc(void *,size_t);														/*	realloc Function but check if the memory allocation has not returned NULL pointer							*/
void 		safefree(void *);																/*	free Function but check if the memory deallocation has returned a NULL pointer								*/
char *		my_strdup(const char *);														/*	Save a copy of a string																						*/
char *		copy_string(char *,int);														/*	Copy a string and return the copied new string																*/
char *		convertToBinary(long int,int);													/*	Convert a Motif Size in a binary equivalent and save it as a character array of the same size 				*/
FASTAFILE *	OpenFasta(char *); 																/* Open a FASTA formatted file for reading																		*/
int			ReadFasta(FASTAFILE *,char **,char **,int *); 									/* Read the sequences from a FASTA formatted file opened														*/
int			WriteFasta(char *,FastaSequences *,int);										/*	Write FASTA sequences from the FataSequences (array of Fasta_seq structures) in file 						*/
void		CloseFasta(FASTAFILE *); 														/* Close a FASTA formatted file opened			 																*/
void		LoadFasta(char *,FastaSequences *,int); 										/* Load the sequences contained in a FASTA formatted file in an array of structure (cf FastaSequences)			*/
char **		create_degeneracy_table(int);													/*	Create the Degeneracy table that will enumerates motifs with wildcards										*/
char **		degenerate_motif(char *,char **,int,int);										/*	Degenerate the undefined position (all positions except extremities) of a motif 							*/
int			verif_X(char *,int);															/*	Check if a Motif contains the character 'X'																	*/
int			compare_motifs(char *,char *);													/*	Compare two Motifs with or without wildcards ( a wildcard match every single letter from the alphabet )		*/
int			motif_scan(Pvoid_t,Pvoid_t *,char *,int,int,char*,FastaSequences *);	/*	Scan the proteome with a list of motifs, and count their occurences in the Proteome							*/
int			parse_enum(char *,char *,int,Pvoid_t *,int);							/*	Read and store the enumeration file in a Judy Array															*/

uint8_t   Index[MAXLINE];		 /* string to insert */

void * safemalloc(size_t s) 
/* malloc with test on the memory allocation*/
{
	void *mptr = malloc(s);
	if(mptr == NULL) {    exit(ENOMEM);    }
	return mptr;
}

void * safecalloc(int taille,size_t s)
/* calloc with test on the memory allocation */
{
	void *mptr = calloc(taille,s);
	if(mptr == NULL) {    exit(ENOMEM);    }
	return mptr;
}

void * saferealloc(void *mptr,size_t s) {
/* realloc with test on the memory allocation*/
	void *newptr = realloc(mptr,s);
	if(newptr == NULL) {    exit(ENOMEM);    }
	return newptr;
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
/************************************PARSER OF FASTA FORMATTED FILE*************************************/
/*******************************************************************************************************/
FASTAFILE * OpenFasta(char *seqfile){

	FASTAFILE *ffp;
	ffp = safemalloc(sizeof(FASTAFILE));

	ffp->fp = fopen(seqfile,"r");              /* Assume seqfile exists & readable!   */
	if (ffp->fp == NULL){
		fprintf(stderr,"Error! Unable to open the FASTA formatted Sequences file \"%s\".\n",seqfile);
		safefree(ffp);
		exit(1);
	}
	
	if ((fgets(ffp->buffer, MAXSEQ, ffp->fp)) == NULL){ safefree(ffp); return NULL; }
	return ffp;
}

int ReadFasta(FASTAFILE *Fin,char **ret_seq, char **ret_name, int *ret_L){

	char *seq = NULL; char * s=NULL; char *name=NULL;
	int nalloc = MAXSEQ;
	int n = 0; 

	seq = safemalloc(sizeof(char) * MAXSEQ);		/* allocate seq in blocks of 10000 residues */
	if (Fin->buffer[0] != '>') return 0;			/* Exit if line buffer does not contain '>' */

	s=strtok(Fin->buffer+1,"\n");
	name=strdup(s);	/* Store the description line */

	while(fgets(Fin->buffer,80,Fin->fp)){
		if (Fin->buffer[0] == '>') break;	/* Next description line was reached */
		Fin->buffer[strlen ( Fin->buffer ) - 1] = '\0';

		for(s=Fin->buffer; *s != '\0'; s++){ /* Store the sequence in an array of char (seq) */
			seq[n]=*s; n++;
			if (nalloc == n){				/* Realloc if there is no room left to store the sequence in the array of char */
				nalloc += MAXSEQ;
				seq = saferealloc(seq, sizeof(char) * nalloc);
			}
		}
	}

	*ret_name = name;
	*ret_seq  = seq;
	*ret_L    = n;
	return 1;
}

int WriteFasta(char * Fastafile, FastaSequences *Sequences, int N){
	int s=0; int pos=0;
	FILE * OUT = NULL;
	
	OUT=fopen(Fastafile,"w");
	for(s=0;s<N;s++){
		fprintf(OUT,">%s\n",Sequences[s]->header);
		for(pos=0;pos<Sequences[s]->Length;pos++){
			if( (pos+1) % 70 == 0 ){ fprintf(OUT,"\n"); }
			fprintf(OUT,"%c",Sequences[s]->sequence[pos]);
		}
		fprintf(OUT,"\n");
	}
	fclose(OUT);
	return 0;
}


void CloseFasta(FASTAFILE *ffp){
	fclose(ffp->fp);
	safefree(ffp);
}

void LoadFasta(char * Fastafile, FastaSequences *Sequences, int N){
/* Takes a path to a fasta file, Loads it into an 2D Array of char */
	int L=0; int Nseq=0;
	char *seq=NULL; char *name=NULL;
	FASTAFILE *IN=NULL;
	
	IN=OpenFasta(Fastafile);
	while(ReadFasta(IN, &seq, &name, &L)){
		Sequences[Nseq]=safemalloc(sizeof(Fasta_seq));
		/*		fprintf(stderr,">%s Length=%d/%d\n", name, L,(int)(strlen(seq)));	*/
		Sequences[Nseq]->header=strdup(name);
/*		fprintf(stderr,"%s\n",  seq);	*/
		Sequences[Nseq]->sequence = safemalloc((L+1)*sizeof(char));
		strncpy(Sequences[Nseq]->sequence,seq,L);
		Sequences[Nseq]->sequence[L]='\0';
		Sequences[Nseq]->Length=(int)(L);
		safefree(seq); safefree(name);
		Nseq++;
	}
	CloseFasta(IN);

}

int verif_X(char * motif,int size){
/* Check if a string contains the character 'X' */
	int i=0;
	for(i=0;i<size;i++){
/*		fprintf(stderr,"%s\n",motif);*/
		if(motif[i] == 'X'){
			return -1;
		}
	}
	return 0;
}

char * copy_string(char * str,int size){
/* Copy a string to a newly allocated character array with user-specified size */
	char *newstr=safemalloc((size+1)*sizeof(char));
	strncpy(newstr,str,size);
	newstr[size] = '\0';
	return newstr;
}

int compare_motifs(char * str1,char * str2){

	int S1=0,S2=0,i=0;
	S1=strlen(str1); S2=strlen(str2);

	if( S1 != S2 ){ 
		return (-1);
	}else{
		for(i=0;i<S1;i++){
/*			fprintf(stderr,"Motif pos %d #%c# =?= Sequence Motif pos %d #%c#\n",i,str1[i],i,str2[i]);*/
			if( str1[i] == str2[i] || str2[i] == '.'){
				continue;
			}else{
				return (1);
			}
		}
		return(0);
	}
}

void my_sprintf(int N,char *string,long int numb){
/* Print a combination (a mask) depending on the size of the motif.*/
	switch (N) {
		case(3):  sprintf(string,"%03ld",numb);  break;
		case(4):  sprintf(string,"%04ld",numb);  break;
		case(5):  sprintf(string,"%05ld",numb);  break;
		case(6):  sprintf(string,"%06ld",numb);  break;
		case(7):  sprintf(string,"%07ld",numb);  break;
		case(8):  sprintf(string,"%08ld",numb);  break;
		case(9):  sprintf(string,"%09ld",numb);  break;
		case(10): sprintf(string,"%010ld",numb); break;
		case(11): sprintf(string,"%011ld",numb); break;
		case(12): sprintf(string,"%012ld",numb); break;
		default:  fprintf(stderr,"Motif size too big (>12 positions)\n"); exit(0);
	}
}

char **degenerate_motif(char *motif,char **combis, int MOTIF_SIZE,int NC) {
/* Degenerate the undefined position (all positions except extremities) of a motif */

	int pos=0; int comb=0; /* int i=0; int cpt=0; int inutile=0; */
	char ** degenerated = safemalloc(NC*sizeof(*degenerated));

/*	fprintf(stderr,"MOTIF TO DEGENERATE: %s | %d combinations\n",motif,NC);*/
	for(comb=0;comb<NC;comb++){
		degenerated[comb] = safemalloc((MOTIF_SIZE+1)*sizeof(**degenerated));
		for(pos=0;pos<MOTIF_SIZE;pos++){
			if(combis[comb][pos] == '0'){
				degenerated[comb][pos]='.';
			} else {
				degenerated[comb][pos]=motif[pos];
			}
		}
		degenerated[comb][MOTIF_SIZE]='\0';
	}
/*	for(i=0;i<NC;i++){	fprintf(stderr,"i=%d\t%s\n",i,degenerated[i]); 	}*/
	return degenerated;
}

char * convertToBinary(long int num,int N) {
/* Convert a Motif Size in a binary equivalent and save it as a character array of the same size */

/*	int i = num;*/
	long int binary = 0; 
	long int place = 0; 
	char *entier= safemalloc((N+1)*sizeof(char));
	entier[N+1] = '\0';
	while (num != 0) { 
		binary = binary + (num%2 * pow(10, place)); 
		num = num / 2; 
		place = place + 1; 
/*		printf("num %d binary %d place %d pow(10, place) %d num mod 2  %d\n",num,binary,place,(int)(pow(10, place)),(num%2));*/
	}
/*	printf("The equivalent of %d binary is %012d\n",i, binary); */
	my_sprintf(N,entier,binary);
	return(entier);
}

char ** create_degeneracy_table(int size){

	int NC=0,i=0;
	char **combis = NULL;					/* 2D character array - Degeneracy table */

	for(i=((pow(2,size)/2)+3);i<pow(2,size);i+=2){
		NC++;
/*		fprintf(stderr,"%d/ size = %d => %s\n",NC,size,convertToBinary((long int)(i),size));*/
	}
	combis = safemalloc( (NC-1) * sizeof(*combis));
/*	 First creates the "degeneracy" table, i.e., it gives all the positions where wildcards are inserted. */
	NC=0;
	for(i=((pow(2,size)/2)+3);i<pow(2,size);i+=2){
		combis[NC]=convertToBinary((long int)(i),size);
		/*printf("i=%d ",i); printf(" %s \t",combis[j]); for(k=0;k<MOTIF_SIZE;k++){	printf(" %c ",combis[j][k]);	}	printf("\n");*/
		NC++;
	}
	return combis;

}
int motif_scan(Pvoid_t k_Array, Pvoid_t * m_Array, char *path_seqfile, int size, int NSeq, char *path_output,FastaSequences * Sequences) {

	int s=0,pos=0,no_X=0,AA1=0,AA2=0,NC=0,mot=0,i=0; /* int cpt=0; */
	char AA[20] = {'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
	char *motif = NULL;						/* character array - Motif */
	char ** combis = NULL;					/* 2D character array - Degeneracy table */
	char ** all_motifs = NULL;				/* 2D character array - Degeneracy table */
	FILE* output = (FILE*)NULL;				/* File Buffer - File Output */
	Pvoid_t   tmpArray = (PWord_t)NULL; 	/* Judy array - Occurences for Motifs within one sequence */
	Pvoid_t   occArray = (PWord_t)NULL; 	/* Judy array - Occurences for Motifs within the whole Proteome*/
	PWord_t   occ,tmp,k,m; 					/* Judy array element */
	Word_t    Rc_word; 						/* size of JudySL array */
/*	int inutile=0;*/

	combis = create_degeneracy_table(size);
	for(i=((pow(2,size)/2)+3);i<pow(2,size);i+=2){
		NC++;
	}

	for(s=0;s<NSeq;s++){
/*		fprintf(stderr,"#seq %d/%d# \n",(s+1),NSeq);*/
		for(pos=0;pos<Sequences[s]->Length-size-1;pos++){
/*			fprintf(stderr,"Seq %d/%d #pos %d/%c (%d)#\n",s+1,NSeq,pos,Sequences[s][pos],(strlen(Sequences[s])-size+1));*/
			motif = copy_string(Sequences[s]->sequence+pos,size);
			no_X=verif_X(motif,size);
			if(no_X == 0){
				all_motifs = degenerate_motif(motif,combis,size,NC);
				for(mot=0;mot<NC;mot++){
					strcpy((char *)(Index),all_motifs[mot]); Index[size]='\0';
					JSLG(k, k_Array,Index);
/*					fprintf(stderr,"%d/%d #mot %s | seq %s# \r",mot,NC,all_motifs[mot],Index);*/
					if(k != NULL) {
/*						fprintf(stderr,"occurences in set = %d \r",(int)(*k));*/
						JSLG(tmp,tmpArray,Index); JSLG(m,*(m_Array),Index); JSLG(occ,occArray,Index);
						if(tmp != NULL ){ /* Motif seen before in the Proteome but within the same sequence */
							(*tmp)++; (*occ)++;
						}else{
							JSLI(tmp, tmpArray,Index);
							(*tmp)++;
							if(occ != NULL && m != NULL){ /* Motif seen before in the Proteome but in a different sequence */
								(*occ)++; (*m)++; 
							}else{						/* Motif never seen before in the Proteome */
								JSLI(m, *(m_Array),Index); JSLI(occ, occArray,Index);
								(*m)++; (*occ)++; 
							}
						}
					}
					safefree(all_motifs[mot]);
				}
				safefree(all_motifs);
			}
/*			fprintf(stderr,"#pos %d : %c - %d \n",pos,Sequences[s][pos],(strlen(Sequences[s])-size+1));*/
			safefree(motif);
		}
		JSLFA(Rc_word, tmpArray);
	}

/*	fprintf(stderr,"%d sequences# \n",(s+1));*/

	output = fopen(path_output, "a");
	if (!output){  fprintf(stderr,"Error! Unable to open the results file \"%s\".\n",path_output); exit(1); }

	for(AA1=0;AA1<20;AA1++) {
		for(AA2=0;AA2<20;AA2++) {
			Index[0] = '\0';
			JSLF(m, *(m_Array), Index);	/* get first string */
			while (m != NULL) {
				JSLG(occ,occArray,Index); 
				if(AA[AA1] == Index[0] && AA[AA2] == Index[size-1]){
					fprintf(output,"%s\t%d\t%d\n", Index,(int)(*occ),(int)(*m));
				}
				JSLN(m, *(m_Array), Index);   /* get next string */
			}
		}
	}
	
	safefree(combis);
	fclose(output);

	JSLFA(Rc_word, k_Array);		/* free array */
	JSLFA(Rc_word, occArray);		/* free array */
	return 0;
}

int parse_enum(char *path_enumfile, char *separateur, int s, Pvoid_t * success_Array, int OnevsAll) {

	FILE* enumfile = (FILE*)NULL;
	char line [MAXLINE]; char *pcourant=NULL;
	PWord_t success;     /* Judy array element. */
/*	int inutile = 0;*/

	enumfile = fopen(path_enumfile, "r");
	/* Check if the file could be open and read */
	if (!enumfile)	{	fprintf(stderr,"Error! Unable to open the CSV formatted file of Enumeration \"%s\".\n",path_enumfile);	exit(1);	}
/*	fprintf(stderr,"%s\n",path_enumfile);*/

	while(fgets(line,MAXHEADER,enumfile)) {
		if( line[0] != '#' ){
/*			fprintf(stderr,"%s\n",line);*/
	/*		scanf("%d",&inutile);*/
			pcourant = strtok(line, separateur);
	/*		printf("pcourant=%s \t",pcourant);*/
			if (pcourant != NULL){
				strcpy((char *)(Index),pcourant); Index[s]='\0';
			}
	/*		printf("Motif=%s \t",(char *)(Index));*/

			pcourant = strtok(NULL, separateur);
			/********************************/
			/* All occurrences per sequence */
			/********************************/
			if(OnevsAll == 0){
				if (pcourant != NULL){
					JSLI(success, *(success_Array),Index);   /* store string into array */
					(*success)=atoi(pcourant);
				}
			}

	/*		printf("pcourant=%s \t",pcourant);*/
			pcourant = strtok(NULL, "\n");
			/*******************************/
			/* One occurrence per sequence */
			/*******************************/
	/*		printf("pcourant=%s \t",pcourant);*/
			if(OnevsAll != 0){
				if (pcourant != NULL){
					JSLI(success, *(success_Array),Index);   /* store string into array */
					(*success)=atoi(pcourant);
				}
			}

/*			fprintf(stderr,"Motif=%s success=%d\n",Index,(int)(*success));*/
	/*		printf("\n");*/
		}
	}
/*	printf("Motif=%s success=%d \t",Index,(int)(*success));*/
	fclose(enumfile);
	return 0;
}

int main (int argc,char **argv){

	/* DECLARATION */
	int NSeq; int size; int i; int OnevsAll;
	char * Prot_seqfile, *Set_enumfile, *Output_path;
	FastaSequences * Sequences;

	Pvoid_t  P_mArray, S_kArray; 
	/* INITIALISATION */

	NSeq=0; size=0; i=0; OnevsAll=0;
	Prot_seqfile = NULL ; Set_enumfile = NULL; Output_path = NULL;
	Sequences=NULL;

	P_mArray = (PWord_t) NULL;  /* Judy array. */
	S_kArray = (PWord_t) NULL;  /* Judy array. */

	/* AFFECTATION */
	Prot_seqfile =  strdup(argv[1]);
	Set_enumfile =  strdup(argv[2]); 
	NSeq         =  atoi(argv[3]); 
	size         =  atoi(argv[4]);
	Output_path  =  strdup(argv[5]);
	if(argc == 7){
		OnevsAll = atoi(argv[6]);
	}
/*	printf("\n");*/
/*	fprintf(stderr,"INPUT :\n Sequences => %s \n Proteome => %s \n",Set_enumfile,Prot_seqfile);*/

	Sequences = safemalloc(NSeq*sizeof(FastaSequences));
	LoadFasta(Prot_seqfile,Sequences,NSeq);
	parse_enum(Set_enumfile,"\t",size,&S_kArray,OnevsAll);
	/* Scan each motif from the set in the Proteome and write a Motif enumeration file for Proteome */
	motif_scan(S_kArray,&P_mArray,Prot_seqfile,size,NSeq,Output_path,Sequences);

	for(i=0;i<NSeq;i++){
		safefree(Sequences[i]->header);
		safefree(Sequences[i]->sequence);
		safefree(Sequences[i]);
	}
	safefree(Sequences);

	safefree(Prot_seqfile);
	safefree(Set_enumfile); 
	safefree(Output_path); 

/*	printf("\n");*/
	return EXIT_SUCCESS;
}

