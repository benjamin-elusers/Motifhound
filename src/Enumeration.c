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
/*                                    Motif Enumeration Process                                  */
/*============================================================================================*/
/* Takes a FASTA file as argument, 4 Integers ( Nb of sequences, Size of motifs, max Nb of Wildcards allowed, Nb of minimal occurences for each motif )  */
/* Loads the FASTA file into Tab of char* , enumerates all motifs and write all of those found more than K_OCC times into a file. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <Judy.h>

#define MIN(a,b)   (a < b ? a : b)
#define MAX(a,b)   (a > b ? a : b)
#define MAXLINE 512 /* max string (line) length */
#define MAXSEQ 40000 /* max sequence length */
#define MAXHEADER 512 /* max header length */
#define JUDYERROR_SAMPLE 1
#define LINELENGTH 50
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

uint8_t   Index[MAXLINE]; /* string to insert */

void *		safemalloc(size_t); 													/*	malloc Function but check if the memory allocation has not returned NULL pointer					*/
void *		safecalloc(int, size_t); 												/*	calloc Function but check if the memory allocation has not returned NULL pointer					*/
void *		saferealloc(void *,size_t);												/*	realloc Function but check if the memory allocation has not returned NULL pointer					*/
void 		safefree(void *);														/*	free Function but check if the memory deallocation has returned a NULL pointer						*/
char *		my_strdup(const char *);												/*	Save a copy of a string																				*/
char *		copy_string(char *,int);												/*	Copy a string and return the copied new string														*/

int			verif_X(char *,int);													/*	Check if a motif contains a 'X'																		*/
void		my_sprintf(int,char *,long int);										/*	Store an integer as a character array and a size equal to N											*/
char *		convertToBinary(long int,int); 											/*	Convert a Motif Size in a binary equivalent and save it as a character array of the same size		*/
FASTAFILE *	OpenFasta(char *);														/*	Open a FASTA formatted file for reading																*/
int			ReadFasta(FASTAFILE *,char **,char **,int *);							/*	Read the sequences from a FASTA formatted file opened												*/
int			WriteFasta(char *,FastaSequences *,int);								/*	Write FASTA sequences from the FataSequences (array of Fasta_seq structures) in file 				*/
void		CloseFasta(FASTAFILE *);												/*	Close a FASTA formatted file opened																	*/
void		LoadFasta(char *,FastaSequences *,int);									/*	Load the sequences contained in a FASTA formatted file in an array of structure (cf FastaSequences)	*/
char **		degenerate_motif(char *,char **,int, int);								/*	Degenerate the undefined position (all positions except extremities) of a motif						*/
int			binomial(int,int);														/*	Calculate binomial coefficient for large numbers													*/
void		Enumeration(int,int,int,char,char,FastaSequences *,char **,int,FILE *);

void * safemalloc(size_t s) {
 /* malloc with test on the memory allocation */
	void *mptr = malloc(s);
	if(mptr == NULL) {    exit(ENOMEM);    }
	return mptr;
}

void * safecalloc(int taille,size_t s) {
/* calloc with test on the memory allocation */
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
	if(p != NULL){
		free(p);
	}
	p=NULL;
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

	if(Nseq != N){ fprintf(stderr,"Input number of sequences not equal to number of fasta sequences in FASTA File.\n"); exit(-2);}
}


void my_sprintf(int N,char *string,long int numb){
/* Print a combination (a mask) depending on the size of the motif.*/
	switch (N) {
		case(3) : sprintf(string,"%03ld",numb);  break;
		case(4) : sprintf(string,"%04ld",numb);  break;
		case(5) : sprintf(string,"%05ld",numb);  break;
		case(6) : sprintf(string,"%06ld",numb);  break;
		case(7) : sprintf(string,"%07ld",numb);  break;
		case(8) : sprintf(string,"%08ld",numb);  break;
		case(9) : sprintf(string,"%09ld",numb);  break;
		case(10): sprintf(string,"%010ld",numb); break;
		case(11): sprintf(string,"%011ld",numb); break;
		case(12): sprintf(string,"%012ld",numb); break;
		case(13): sprintf(string,"%013ld",numb); break;
		case(14): sprintf(string,"%014ld",numb); break;
		case(15): sprintf(string,"%015ld",numb); break;
		case(16): sprintf(string,"%016ld",numb); break;
		case(17): sprintf(string,"%017ld",numb); break;
		case(18): sprintf(string,"%018ld",numb); break;
		case(19): sprintf(string,"%019ld",numb); break;
		case(20): sprintf(string,"%020ld",numb); break;
		default : fprintf(stderr,"Problem with the combinations size\n"); exit(-N);
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
	entier[N] = '\0';
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

char * copy_string(char * str,int size){
/* Copy a string to a newly allocated character array with user-specified size */
	char *newstr=safemalloc((size+1)*sizeof(char));
	strncpy(newstr,str,size);
	newstr[size] = '\0';
	return newstr;
}

int verif_X(char * motif,int size){
	
	int i=0;
	for(i=0;i<size;i++){
		if(motif[i] == 'X'){
			return -1;
		}
	}
	return 0;
}

void Enumeration (int MOTIF_SIZE,int N,int K_OCC,char AA1,char AA2,FastaSequences * Sequences,char **combis, int NC,  FILE * OUT){

	int k=0; int pos=0; int s=0; int no_X=-1; /* int p=0; */

	char* motif = NULL;
	char ** all_motifs = NULL; 

	Pvoid_t   occArray = (PWord_t)NULL;  /* Judy array. */
	Pvoid_t   seqArray = (PWord_t)NULL;  /* Judy array. */
	Pvoid_t   tmpArray = (PWord_t)NULL;  /* Judy array. */

	PWord_t   occ;		     /* Judy array element. */
	PWord_t   seq;		     /* Judy array element. */
	PWord_t   tmp;		     /* Judy array element. */

	Word_t    Rc_word;		     /* size of JudySL array. */
	for(s=0;s<N;s++){
/*		printf("SEQ #%d# = %s\n",(s+1),Sequences[s]);*/
		for(pos=0;pos<Sequences[s]->Length-MOTIF_SIZE-1;pos++){
			if(Sequences[s]->sequence[pos] == AA1 && Sequences[s]->sequence[pos+MOTIF_SIZE-1] ==  AA2){
				motif = copy_string(Sequences[s]->sequence+pos,MOTIF_SIZE);
				no_X=verif_X(motif,MOTIF_SIZE);
				if(no_X == 0){
/*					fprintf(stderr," Motif %s + First (Seq) %c == (fixed AA) %c && Last (Seq) %c == (fixed AA) %c +\n",motif,Sequences[s]->sequence[pos],AA1,Sequences[s]->sequence[pos+MOTIF_SIZE],AA2);*/
/*					fprintf(stderr," Motif %s pos %d %c\n",motif,pos,Sequences[s]->sequence[pos]);*/
					all_motifs = degenerate_motif(motif,combis,MOTIF_SIZE,NC);
					for(k=0;k<NC;k++){
	/*					fprintf(stderr,"k=%d\t%s\n",k,all_motifs);*/
						strcpy((char *)(Index),all_motifs[k]); Index[MOTIF_SIZE]='\0';
						JSLG(tmp,tmpArray,Index); JSLG(occ,occArray,Index); JSLG(seq,seqArray,Index);
						if(tmp != NULL ){
							(*tmp)++; (*occ)++;
						}else{
							if(occ != NULL && seq != NULL){
								(*occ)++; (*seq)++; 
							}else{
								JSLI(seq, seqArray, Index);   /* store string into array */
								(*seq)++;
								JSLI(occ, occArray,Index);   /* store string into array */
								(*occ)++; 
							}
							JSLI(tmp, tmpArray,Index);   /* store string into array */
							(*tmp)++;
						}
						safefree(all_motifs[k]); 		/* free array */
					}
					safefree(motif);
					safefree(all_motifs); 		/* free array */
				}
			}
		}
		JSLFA(Rc_word, tmpArray);
	}
/*	fprintf(stderr,"\n");*/

	/* Prints the motifs sorted lexicographically */
	Index[0] = '\0';
	JSLF(seq, seqArray, Index);	/* get first string */
	while (seq != NULL) {
		JSLG(occ, occArray, Index);	/* get first string */
		if(occ != NULL) {
			if((int)(*seq) >= K_OCC){
				fprintf(OUT,"%s\t%d\t%d\n", Index,(int)(*occ),(int)(*seq));
			}
		}
		JSLN(seq, seqArray, Index);   /* get next string */
	}
	JSLFA(Rc_word, seqArray);		/* free array */
	JSLFA(Rc_word, occArray);		/* free array */
}

/*Binomial coefficient for large numbers*/
int binomial(int n,int k){
	int r = 1;
	int d = n - k;
	/* choose the smaller of k and n - k */
	if (d > k) { k = d; d = n - k; }

	while (n > k) {
		r *= n--;
		/* divide (n - k)! as soon as we can to delay overflows */
		while (d > 1 && !(r % d)) r /= d--;
	}
	return r;
}

int main(int argc,char **argv){ 

	/* DECLARATIONS & INITIALISATIONS */
	int i=0; int j=0; int k=0; int x=0;
	int AA1=0; int AA2=0;

	int N = 0; int K_OCC = -1;
	int MOTIF_SIZE = -1; int MIN_DEF = -1; int MAX_DEF = -1; int N_DEF=0;
	int MAX_COMB=0; int MIN_COMB=0; int NB_COMB=0; int NB_COMB_DEF=0;

	char *Fastafile = NULL; char *Outdir = NULL; char path_output[256];
	/* Amino acids ordered by Alphabetical order */
	char AA[20] = {'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'}; 
	/* Amino Acids orderered by their Frequency in Yeast */
/*	char AA[20] = {'L','S','K','E','I','V','A','D','T','N','G','R','F','P','Q','Y','H','M','C','W'};*/

	char **combis = NULL; char **masks = NULL;

	FastaSequences *Sequences = NULL;

	FILE * OUT = (FILE*) NULL;

	/* AFFECTATIONS */
	if( argc < 6 ){
		fprintf(stderr,"NOT ENOUGH ARGUMENTS !                                           \n");
		fprintf(stderr,"-----------------------------------------------------------------\n");
		fprintf(stderr,"USAGE  : %s ARG1  ARG2  ARG3  ARG4  ARG5  [ARG6] or [ARG6  ARG7] \n",argv[0]);
		fprintf(stderr,"ARG1   : Location of the Sequence file                           \n");
		fprintf(stderr,"ARG2   : Result file path                                        \n");
		fprintf(stderr,"ARG3   : Number of Sequences                                     \n");
		fprintf(stderr,"ARG4   : Motif Length                                            \n");
		fprintf(stderr,"ARG5   : Minimal number of occurences                            \n");
		fprintf(stderr,"-----------------------------------------------------------------\n");
		fprintf(stderr,"[ARGx] : Optional Arguments                                      \n");
		fprintf(stderr,"ARG6   : MINIMUM number of non-wildcard positions                \n");
		fprintf(stderr,"ARG7   : MAXIMUM number of non-wildcard positions                \n");
		exit(-argc);
	}

	Fastafile = strdup(argv[1]);
	Outdir = strdup(argv[2]);
	N = atoi(argv[3]);
	MOTIF_SIZE = atoi(argv[4]);
	K_OCC = atoi(argv[5]);


	if( argc > 6 ){
		if( argc==7 ){ MIN_DEF = atoi(argv[6]); MAX_DEF = MOTIF_SIZE;    }
		if( argc==8 ){ MIN_DEF = atoi(argv[6]); MAX_DEF = atoi(argv[7]); }
	}else{
		MIN_DEF=3;
		MAX_DEF=MOTIF_SIZE;
	}
	
	/* PRINT OPTIONS */
	sprintf(path_output,"%s",Outdir);
/*	printf("\n\nFASTA Input = %s\nEnum Output = %s\n\n",Fastafile,path_output);*/
/*	printf("Nseq = %d   ",N);*/
/*	printf("MOTIF_SIZE =%d   ",MOTIF_SIZE);*/
/*	printf("K_OCC =%d\t\n\n",K_OCC);*/
/*	printf("N_WILD = %d\t\n",N_WILD);*/
	if(K_OCC != -1 && MOTIF_SIZE != -1 && MIN_DEF != -1 && MAX_DEF != -1){
		MIN_COMB= pow(2,MOTIF_SIZE-1) + pow(2,0); MAX_COMB= pow(2,MOTIF_SIZE-1) + pow(2,0);
		for(i=1;i<=1;i++){ MIN_COMB += pow(2,i); }
		for(i=1;i<MOTIF_SIZE-1;i++){ MAX_COMB += pow(2,i); }
		for(i=MIN_COMB;i<=MAX_COMB;i+=2){ NB_COMB += 1; }
		for(i=MIN_DEF-2;i<=MAX_DEF-2;i++){ NB_COMB_DEF += binomial(MOTIF_SIZE-2,i); }
/*		printf("MOTIF_SIZE = %d\t",MOTIF_SIZE);*/
/*		printf("MIN_DEF = %d\t",MIN_DEF);*/
/*		printf("MAX_DEF = %d\n",MAX_DEF);*/
/*		printf("MIN_COMB = %d (%s)\t",MIN_COMB,convertToBinary((long int)(MIN_COMB),MOTIF_SIZE));*/
/*		printf("MAX_COMB = %d (%s)\t",MAX_COMB,convertToBinary((long int)(MAX_COMB),MOTIF_SIZE));*/
/*		printf("NB_COMB = %d\t\n",NB_COMB);*/
/*		printf("NB_COMB_DEF = %d\t\n",NB_COMB_DEF);*/

/*		printf("b(5,1) = %d\t\n",binomial(5,1));*/
/*		printf("b(5,2) = %d\t\n",binomial(5,2));*/
/*		printf("b(5,3) = %d\t\n",binomial(5,3));*/
/*		printf("b(5,4) = %d\t\n",binomial(5,4));*/
/*		printf("b(5,5) = %d\t\n",binomial(5,5));*/

		combis = safemalloc(NB_COMB*sizeof(*combis));
		masks = safemalloc(NB_COMB_DEF*sizeof(*masks));
		Sequences = safemalloc(N*sizeof(FastaSequences));
		LoadFasta(Fastafile,Sequences,N);
/*		WriteFasta("/home/benjamin/Desktop/Protein_Short_Linear_Motifs/Articles/Functional_Datasets/Construction_of_human_activity_based_phosphorylation_networks/Results/Enumeration/MAPK1/Data_User_Input/Test.faa",Sequences,N);*/
/*		for(i=0;i<N;i++){ 	fprintf(stderr,"ind=%d len %d seq=%s\n",i,(int)(strlen(Sequences[i])),Sequences[i]);	} */

/*		 First creates the "degeneracy" table, i.e., it gives all the positions where wildcards are inserted. */
		i=0;
		for(j=MIN_COMB;j<=MAX_COMB;j+=2){
			N_DEF=0;
			combis[i]=convertToBinary((long int)(j),MOTIF_SIZE);
			for(k=0; k<MOTIF_SIZE;k++){ N_DEF+=combis[i][k]-'0'; } /* Substract the ascii value of char '0' */
/*			printf("j=%d ",j); for(k=0;k<MOTIF_SIZE;k++){	printf(" %c ",comb[k]);	}	printf("\nN_DEF=%d\n",N_DEF); */
			if( N_DEF >= MIN_DEF && N_DEF <= MAX_DEF ){ masks[x]=combis[i]; x++; }
			i++;
		}
/*		printf("Number of combinations : %d\n",NB_COMB);*/
/*		printf("Number of combinations ( %d <= # non-wildcard positions <= %d ) : %d\n",MIN_DEF,MAX_DEF,NB_COMB_DEF);*/
		safefree(combis);

		OUT=fopen(path_output,"a");
/*		fprintf(stderr,"#");*/
		for(AA1=0;AA1<20;AA1++) {
/*			fprintf(stderr," -%c- ",AA[AA1]);*/
			for(AA2=0;AA2<20;AA2++) {
				Enumeration(MOTIF_SIZE,N,K_OCC,AA[AA1],AA[AA2],Sequences,masks,NB_COMB_DEF,OUT);
			}
/*			fprintf(stderr," %c ",AA[AA1]);	*/
		}
/*		fprintf(stderr,"\n");	*/
		fclose(OUT);
		
		for(i=0;i<N;i++){
			safefree(Sequences[i]->header);
			safefree(Sequences[i]->sequence);
			safefree(Sequences[i]);
		}
		safefree(Sequences);
		safefree(masks);
		safefree(Fastafile);
		safefree(Outdir);

/*		fprintf(stderr,"Enumeration Done#\n");	*/
	}
	return(0);
}

