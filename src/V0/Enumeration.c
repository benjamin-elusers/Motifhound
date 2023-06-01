/* +------------------------------------------------------------+ */
/* |  Date		    : 16/05/2011                                  | */
/* |  Last Update : 28/05/2023                                  | */
/* |  Author      : B. Dubreuil                                 | */
/* |  Contact		  : dubreuil.benjamin@hotmail.com               | */
/* |  Affiliation : Michnik lab - Université de Montréal (UDEM) | */
/* |  Centre de Bioinformatique et Genomique Robert Cedergren   | */
/* +------------------------------------------------------------+ */
/*============================================================================================*/
/*##     ##  #######  ######## #### ######## ##     ##  #######  ##     ## ##    ## ########  */
/*###   ### ##     ##    ##     ##  ##       ##     ## ##     ## ##     ## ###   ## ##     ## */
/*#### #### ##     ##    ##     ##  ##       ##     ## ##     ## ##     ## ####  ## ##     ## */
/*## ### ## ##     ##    ##     ##  ######   ######### ##     ## ##     ## ## ## ## ##     ## */
/*##     ## ##     ##    ##     ##  ##       ##     ## ##     ## ##     ## ##  #### ##     ## */
/*##     ## ##     ##    ##     ##  ##       ##     ## ##     ## ##     ## ##   ### ##     ## */
/*##     ##  #######     ##    #### ##       ##     ##  #######   #######  ##    ## ########  */
/*                                     Motif Enumeration                                      */
/*============================================================================================*/
/* Loads FASTA fromatted sequences from a file */
/* Enumerate sequence motifs with user input on (i) specific length (ii) number of wildcards (iii) minimum occurrence per motif */

/* Compile the following script */
/* gcc  Enumeration.c safealloc.c fasta.c motif.c -o Enumeration.exe -lm -lJudy -lpthread */ 

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h> 
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "Judy.h"
#include "safealloc.h"
#include "motif.h"
#include "fasta.h"

#define MIN(a,b)   (a < b ? a : b)
#define MAX(a,b)   (a > b ? a : b)
#define JUDYERROR_SAMPLE 1
#define LINELENGTH 50

uint8_t   Index[MAXLINE]; /* string to insert */

typedef struct {
    char* fastafile;
    char* output;
    char* seqtype;
    int mlen;
    int kmin;
    int xmin;
    int xmax;
} OptionValues;


void printUsage (); /* Print the usage of the program listing all the options */
int  checkOptionValues (const struct option*, OptionValues *); /* Validate option values */
int  sumArray (const int, int); /*	Calculate sum of values in an integer array */
int	 binomial (int, int);	/*	Calculate binomial coefficient for large numbers */
void Enumeration (int, char **, int, int, char, char, SEQUENCES *, int,FILE *);


void printUsage() {
	printf("Usage: Enumeration.exe -f <fastafile> -o <output> -m <mlen> -k <kmin> [-x1 <xmin>] [-x2 <xmax>] [-a <seqtype>]\n");
	printf("Arguments (* = required):\n");
	printf("* -f, --fasta     <fastafile>   Path to the input FASTA file\n");
	printf("* -o, --output    <outputfile>  Path to the output file\n");
	printf("* -a, --seqtype   <seqtype>     Sequence type (AA or NT)\n");
	printf("* -m, --mlen      <mlen>        Motif length\n");
	printf("  -k, --kmin      <kmin>        Minimal number of occurrences\n");
	printf("  -x, --xmin      <xmin>        Minimum number of non-wildcard positions\n");
	printf("  -y, --xmax      <xmax>        Maximum number of non-wildcard positions\n");
}


int checkOptionValues(const struct option* long_options,OptionValues *options_values) {

	for (int i = 0; long_options[i].name != NULL; i++) {
		/*printf("Option: %s (short %c)\n", long_options[i].name,long_options[i].val); */
		/*printf("required: %s\n", long_options[i].has_arg == required_argument ? "Yes" : "No"); */

		if(options_values->fastafile == NULL){
			fprintf(stderr,"Missing fastafile!\n");
			exit(404);
		}
		
		if(options_values->output == NULL){
			fprintf(stderr,"Missing output file!\n");
			exit(404);
	}
		
		if(options_values->seqtype == NULL){
			fprintf(stderr,"Missing type of fasta sequences!\n");
			exit(404);
		}
		
		if(options_values->mlen == -1){
			fprintf(stderr,"Missing motif length!\n");
			exit(404);
		}

		if(options_values->kmin == -1){
			fprintf(stderr,"Missing minimum of motif occurrences (default to 1)!\n");
			exit(404);
		}

		if(options_values->xmin == -1 || options_values->xmin < 3){
			fprintf(stderr,"Missing minimum non-wildcard positions (default to 3)!\n");
			options_values->xmin = 3;
		}

		if(options_values->xmax == -1 || options_values->xmax > options_values->mlen){
			fprintf(stderr,"Missing maximum non-wildcard positions (default to motif size)!\n");
			options_values->xmax = options_values->mlen;
		}
		
	}

	fprintf(stderr,"fastafile is set to : %s\n",options_values->fastafile);
	fprintf(stderr,"output is set to    : %s\n",options_values->output);
	fprintf(stderr,"seqtype is set to   : %s\n",options_values->seqtype);
	fprintf(stderr,"mlen is set to      : %d\n",options_values->mlen);
	fprintf(stderr,"kmin is set to      : %d\n",options_values->kmin);
	fprintf(stderr,"xmin is set to      : %d\n",options_values->xmin);
	fprintf(stderr,"xmax is set to      : %d\n",options_values->xmax);
	
	
	/* Check if required options are provided */
	if (options_values->fastafile == NULL || options_values->seqtype == NULL || options_values->output == NULL || options_values->mlen == -1) {
		printUsage();
		exit(404);
	}
	return 0;
}


/*Sum integer values in an array*/
int sumArray(const int *array, int length) {
    int sum = 0;
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }
    return sum;
}


/*Binomial coefficient for large numbers*/
int binomial(int n,int k){
	int r = 1;
	int d = n - k;
	/* choose the smaller of k and n - k */
	if (d > k) {
		k = d; 
		d = n - k;
	}

	while (n > k) {
		r *= n--;
		/* divide (n - k)! as soon as we can to delay overflows */
		while (d > 1 && !(r % d)){
			r /= d--;
		}
	}
	return r;
}




void Enumeration (int len, char **masks, int NC, int K_OCC, char firstletter, char lastletter, SEQUENCES * sequences, int nseq,  FILE * OUT){

	int k=0; int pos=0; int s=0; 

	char* motif = NULL;
	char ** all_motifs = NULL; 

	Pvoid_t   occArray = (PWord_t)NULL;  /* Judy array. */
	Pvoid_t   seqArray = (PWord_t)NULL;  /* Judy array. */
	Pvoid_t   tmpArray = (PWord_t)NULL;  /* Judy array. */

	PWord_t   occ;				 /* Judy array element. */
	PWord_t   seq;				 /* Judy array element. */
	PWord_t   tmp;				 /* Judy array element. */

	Word_t		Rc_word;				 /* size of JudySL array. */
	for(s=0;s<nseq;s++){
		for(pos=0;pos<sequences[s]->length-len-1;pos++){
			if(sequences[s]->sequence[pos] == firstletter && sequences[s]->sequence[pos+len-1] ==  lastletter){
				motif = CopyString(sequences[s]->sequence+pos,len);
				if(VerifyX(motif,len) == 0){
					all_motifs = DegenerateMotif(masks,motif,len,NC);
					for(k=0;k<NC;k++){
						strcpy((char *)(Index),all_motifs[k]); Index[len]='\0';
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
						safefree((void**)&all_motifs[k]); 		/* free array */
					}
					safefree((void**)&motif);
					safefree((void**)&all_motifs); 		/* free array */
				}
			}
		}
		JSLFA(Rc_word, tmpArray);
	}

	/* Prints the motifs sorted lexicographically */
	Index[0] = '\0';
	JSLF(seq, seqArray, Index);	/* get first string */
	while (seq != NULL) {
		JSLG(occ, occArray, Index);	/* get first string */
		if(occ != NULL) {
			if((int)(*seq) >= K_OCC){
				/*fprintf(stderr,"%s\t%d\t%d\n", Index,(int)(*occ),(int)(*seq));*/
				fprintf(OUT,"%s\t%d\t%d\n", Index,(int)(*occ),(int)(*seq));
			}
		}
		JSLN(seq, seqArray, Index);   /* get next string */
	}
	JSLFA(Rc_word, seqArray);		/* free array */
	JSLFA(Rc_word, occArray);		/* free array */
}


int main(int argc,char **argv){ 

	/* DECLARATIONS & INITIALISATIONS */
	int i=0; int j=0; int k=0; int x=0; int nx=0;
	int first=0; int last=0;
	int nseq = 0; int cmax=0; int cmin=0; int nc=0; int ncx=0;
	int opt; int nletters=0;
	int * cdef = NULL;

	char *letters = NULL; char *outdir = NULL;
	char **combinations = NULL; char **masks = NULL;
	SEQUENCES *sequences = NULL;
	FILE * OUT = (FILE*) NULL;

	/* Variables to store the option values */
	OptionValues options = {NULL, NULL, NULL, -1, 3, 3, 3};
	
	const char *short_options = "f:o:m:k:x:y:a:";
	/* Define the long options array */
	const struct option long_options[] = {
		{"fasta", required_argument, NULL, 'f'},
		{"output", required_argument, NULL, 'o'},
		{"seqtype", required_argument, NULL, 'a'},
		{"mlen", required_argument, NULL, 'm'},
		{"kmin", optional_argument, NULL, 'k'},
		{"xmin", optional_argument, NULL, 'x'},
		{"xmax", optional_argument, NULL, 'y'},
		{NULL, 0, NULL, 0}
	};
		

	options.xmin=3;
	options.kmin=1;
	printf("-------------------INPUT-------------------\n");
	while ((opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
		switch (opt) {
			case 'f':
				options.fastafile = strdup(optarg);
				break;
			case 'o':
				options.output = strdup(optarg);
				break;
			case 'a':
					options.seqtype = optarg;
					break;
			case 'm':
					options.mlen = atoi(optarg);
					options.xmax=options.mlen;
					break;
			case 'k':
					if (optarg != NULL) {
						options.kmin = atoi(optarg);
					}
					break;
			case 'x':
					if (optarg != NULL) {
							options.xmin = atoi(optarg);
					}
					break;
			case 'y':
					if (optarg != NULL) {
							options.xmax = atoi(optarg);
					}
					break;
			case '?':
					printf("Unknown option: %c\n", optopt);
					break;
			default:
				printUsage();
				return 1;		
		}
	}
	checkOptionValues(long_options,&options);
	printf("---------------------------------------------\n");


	// Process non-option arguments (if any)
	for (i = optind; i < argc; i++) {
		printf("Non-option argument: %s\n", argv[i]);
	}

	if( strcmp(options.seqtype, "AA") == 0 ){
		letters = strdup(alphabet.AA);
	}else if( strcmp(options.seqtype, "NT") == 0 ){
		letters = strdup(alphabet.NT);
	}

	nletters = strlen(letters);
	fprintf(stderr,"sequence alphabet has %d letters : [%s]\n",nletters,letters);

	cmin = pow(2,options.mlen-1) + pow(2,1) + pow(2,0); 
	cmax = cmin;
	for(i=2;i<options.mlen-1;i++){ cmax += pow(2,i); }
	fprintf(stderr,"cmin %d \n",cmin);
	fprintf(stderr,"cmax %d \n",cmax);
	fprintf(stderr,"# Combinations %d \n",cmax-cmin);

	nc = ((cmax-cmin)/2) + 1;
	fprintf(stderr,"# Masks %d \n",nc);

	cdef = safecalloc(0,options.mlen);
	fprintf(stderr,"binomial\n");
	for(i=options.xmin-2;i<=options.xmax-2 & i <=options.mlen-2;i++){
		cdef[i] = binomial(options.mlen-2,i)  ;
		/* fprintf(stderr,"%d=%d\n",i+2,cdef[i]); */
	}
	ncx = sumArray(cdef,options.mlen);
	fprintf(stderr,"# selected Masks %d \n",ncx);


	combinations = safemalloc(nc*sizeof(*combinations));
	masks = safemalloc(ncx*sizeof(*masks));
	/* First creates the "degeneracy" table, i.e., it gives all the positions where wildcards are inserted. */
	i=0;
	for(j=cmin;j<=cmax;j+=2){
		nx=0;
		combinations[i]=ToBinary((unsigned long long int)(j),options.mlen);
		/*fprintf(stderr,"%02d) %s\n",i+1,combinations[i]);*/
		for(k=0; k<options.mlen;k++){ nx +=combinations[i][k]-'0'; } /* Substract the ascii value of char '0' */
		/*fprintf(stderr,"number of 1s = %d\n",nx);*/
		if( nx >= options.xmin && nx <= options.xmax ){ 
			masks[x]=combinations[i]; 
			fprintf(stderr,"VALID MASK %d. %s\n",x+1,masks[x]);
			x++; 
		}
		i++;
	}
	safefree((void**)&combinations);
	fprintf(stderr,"Number of degenerate combinations is: %d\n",nc);
	fprintf(stderr,"Wildcard combination from %s to %s\n",ToBinary(cmin,options.mlen),ToBinary(cmax,options.mlen));
	fprintf(stderr,"Number of degenerate combinations between %d and %d non-wildcard positions: %d\n",options.xmin,options.xmax,ncx);

	nseq = countSequences(options.fastafile);
	fprintf(stderr,"Number of sequences found: %d\n",nseq);
	sequences = safemalloc(nseq*sizeof(SEQUENCES));
	LoadFasta(options.fastafile,sequences,nseq);


	fprintf(stderr, "output = %s\n",options.output);
	OUT=fopen(options.output,"a");
	if (OUT == NULL) {
		printf("Error opening output file.\n");
		return -1;
	}

	for(first=0;first<nletters;first++) {
		for(last=0;last<nletters;last++) {
			fprintf(stderr,"first %c last %c\n", letters[first], letters[last]);
			Enumeration(options.mlen, masks, ncx, options.kmin, letters[first], letters[last], sequences, nseq, OUT);
		}
	}
	fclose(OUT);
	
	fprintf(stderr,"done Enumeration\n");

	for(i=0;i<nseq;i++){
		safefree((void**)&sequences[i]->header);
		safefree((void**)&sequences[i]->sequence);
		safefree((void**)&sequences[i]);
	}
	safefree((void**)&sequences);
	safefree((void**)&masks);
	safefree((void**)&options.fastafile);
	safefree((void**)&outdir);

	return(0);
}