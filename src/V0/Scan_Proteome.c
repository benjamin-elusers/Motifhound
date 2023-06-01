/* +------------------------------------------------------------+ */
/* |  Date        : 10/05/2012                                  | */
/* |  Last Update : 14/05/2023                                  | */
/* |  Author      : B. Dubreuil                                 | */
/* |  Contact     : dubreuil.benjamin@hotmail.com               | */
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
#define MOTIFMAXLEN 30
#define FASTALINEWIDTH 80 /* max number of symbols in fasta formatted sequence */
#define HEADERMAXLEN 512 /* max header length */
#define SEQMAXLEN 40000 /* max sequence length */
#define MAXLINE 512 /* max string length */
#define strdup(x) DuplicateString(x)

const char AMINOACIDS[] = "ACDEFGHIKLMNPQRSTUVWY";
const char BASES[] = "ACGTU";

typedef struct {
  FILE* fp;
  char  buffer[SEQMAXLEN];
} FASTAFILE;

typedef struct {
	char* header;
	char* sequence;
	int length;
} SEQUENCE;

typedef SEQUENCE* SEQUENCES;

uint8_t   Index[MAXLINE];		 /* string to insert */

/* Memory allocation with test allocation (not a NULL pointer) */
void*  safemalloc(size_t);         // Allocate memory without initialization of values
void*  safecalloc(int, size_t);    // Allocate and initialize memory to 0
void*  saferealloc(void*,size_t);  // Reallocate memory
void   safefree(void*);            // Deallocate memory safely (set to null pointer)

/* Functions for string manipulation */
char*  DuplicateString(const char*); // Copy a string identically
char*  CopyString(char*,int);        // Copy the first n characters of a string
char*  ToBinary(long int,int);	     // Convert long integer into a binary format

/* Functions to manipulate FASTA formatted sequences */
FASTAFILE*	OpenFasta(char*);                    // Open FASTA file for reading
int    ReadFasta(FASTAFILE*,char**,char**,int*); // Read sequences from FASTA file
int    WriteFasta(char*,SEQUENCES*,int);         //	Write sequences to output FASTA file
void   CloseFasta(FASTAFILE*);                   // Close FASTA file
void   LoadFasta(char*,SEQUENCES*,int);          // Load sequences from FASTA formatted file into an array of structure (SEQUENCES)

/* Functions for motifs */
char** MaskMotifPositions(int);               // Create masks of degenerate motifs (i.e. with wildcards positions)
char** DegenerateMotif(char*,char**,int,int); // Degenerate motifs positions according to the masks
int    VerifyX(char*,int);                    // Check whether a motif contains the character 'X'
int    CompareMotifs(char *,char *);          // Compare two strings

/* Main functions */
int    ScanMotif(Pvoid_t,Pvoid_t*,char*,int,int,char*,SEQUENCES*);	/*	Scan the proteome with a list of motifs, and count their occurences in the Proteome							*/
int    ReadEnumeration(char*,char*,int,Pvoid_t *,int);			/*	Read and store the enumeration file in a Judy Array															*/


void* safemalloc(size_t num, size_t size) {
/* calloc (initialization to 0) with test on memory allocation */
    void* ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(ENOMEM);
    }
    return ptr;
}

void* safecalloc(size_t num, size_t size) {
/* calloc (initialization to 0) with test on memory allocation */
    void* ptr = calloc(num, size);
    if (ptr == NULL) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(ENOMEM);
    }
    return ptr;
}

void* saferealloc(void *ptr,size_t s) {
/* realloc with test on successful memory allocation */
	void *newptr = realloc(ptr,s);
	if(newptr == NULL) { 
		fprintf(stderr, "Memory reallocation failed!\n");
		exit(ENOMEM);    }
	return newptr;
}

void safefree(void** ptr){
     free(*ptr);
    *ptr = NULL;
}


/**
 * Open a FASTA formatted file.
 *
 * This function opens a FASTA formatted file for reading. It returns a pointer
 * to a FASTAFILE structure on success, and NULL on failure.
 *
 * @param seqfile Path to the FASTA formatted file.
 * @return Pointer to a FASTAFILE structure on success, NULL on failure.
 */
FASTAFILE* OpenFasta(const char* seqfile){
    // Allocate memory for the FASTAFILE structure
    FASTAFILE* ffp = safemalloc(sizeof(FASTAFILE));

    // Try to open the file
    ffp->fp = fopen(seqfile, "r");

    // If the file could not be opened, print an error message, free the memory
    // allocated for the FASTAFILE structure, and return NULL
    if (ffp->fp == NULL){
        fprintf(stderr, "Error! Unable to open the FASTA formatted Sequences file \"%s\".\n", seqfile);
        safefree((void**)&ffp);
        return NULL;
    }

    // Try to read the first line of the file
    // If the line could not be read, free the memory allocated for the FASTAFILE
    // structure and return NULL
    if ((fgets(ffp->buffer, SEQMAXLEN, ffp->fp)) == NULL){
        safefree((void**)&ffp);
        return NULL;
    }

    return ffp;
}

/**
 * Read a sequence from a FASTA formatted file.
 *
 * This function reads a single sequence from a FASTA formatted file. It returns 1
 * on success and 0 on failure. On success, it sets the values pointed to by SEQ,
 * DESC, and L to the sequence, the name of the sequence, and the length of
 * the sequence, respectively.
 *
 * @param INPUT A pointer to a FASTAFILE structure.
 * @param SEQ A pointer to a char pointer. On success, *SEQ is set to the sequence.
 * @param DESC A pointer to a char pointer. On success, *DESC is set to the name of the sequence.
 * @param L A pointer to an int. On success, *L is set to the length of the sequence.
 * @return 1 on success, 0 on failure.
 */
int ReadFasta(FASTAFILE* INPUT, char** SEQ, char** DESC, int* L){
    // Initializations
    char* sequence = safemalloc(sizeof(char) * SEQMAXLEN);  // Allocate string in blocks of SEQMAXLEN residues
    char* s = NULL; 
    char* description = NULL;
    int nalloc = SEQMAXLEN;
    int n = 0; 

    // Exit if line buffer does not contain '>'
    if (INPUT->buffer[0] != '>') {
        safefree((void**)&sequence);
        return 0; 
    }

    // Store the description line
    s = strtok(INPUT->buffer + 1, "\n");
    description = strdup(s);	

    // Read the sequence
    while(fgets(INPUT->buffer, FASTALINEWIDTH, INPUT->fp)){
        // If the next description line was reached, stop reading
        if (INPUT->buffer[0] == '>') break;
        
        // Remove the newline character from the end of the buffer
        INPUT->buffer[strlen(INPUT->buffer) - 1] = '\0';

        // Store the sequence in an array of char (sequence)
        for(s = INPUT->buffer; *s != '\0'; s++){ 
            sequence[n] = *s;
            n++;
            
            // If there is no room left to store the sequence in the array of char, realloc
            if (nalloc == n){
                nalloc += SEQMAXLEN;
                sequence = saferealloc(sequence, sizeof(char) * nalloc);
            }
        }
    }

    // Set return values
    *DESC = description;
    *SEQ = sequence;
    *L = n;

    return 1;
}

/**
 * Write sequences to a FASTA formatted file.
 *
 * This function writes sequences to a FASTA formatted file. Each sequence is written
 * with a line width of 70 characters.
 *
 * @param fastafile The path of the FASTA file to write to.
 * @param sequences A pointer to an array of SEQUENCES structures.
 * @param N The number of sequences to write.
 * @return 0 on success, non-zero on failure.
 */
int WriteFasta(char * fastafile, SEQUENCES *sequences, int N) {
    int s = 0;
    int pos = 0;
    FILE * OUT = fopen(fastafile, "w");

    if (OUT == NULL) {
        fprintf(stderr, "Error opening output file %s\n", fastafile);
        return 1;
    }

    for(s = 0; s < N; s++) {
        fprintf(OUT, ">%s\n", sequences[s]->header);
        for(pos = 0; pos < sequences[s]->length; pos++) {
            if ((pos + 1) % FASTALINEWIDTH == 0) {
                fprintf(OUT, "\n");
            }
            fprintf(OUT, "%c", sequences[s]->sequence[pos]);
        }
        fprintf(OUT, "\n");
    }

    fclose(OUT);
    return 0;
}

/**
 * Close a FASTA file and free the associated structure.
 *
 * This function closes a FASTA file and frees the associated FASTAFILE structure.
 *
 * @param ffp A pointer to a FASTAFILE structure.
 */
void CloseFasta(FASTAFILE *ffp) {
    fclose(ffp->fp);
    safefree((void**)&ffp);
}

/**
 * Load sequences from a FASTA formatted file.
 *
 * This function reads sequences from a FASTA formatted file and stores them in an
 * array of SEQUENCES structures.
 *
 * @param fastafile The path of the FASTA file to read from.
 * @param sequences A pointer to an array of SEQUENCES structures.
 * @param N The maximum number of sequences to read.
 */
void LoadFasta(char * fastafile, SEQUENCES *sequences, int N) {
    int L = 0, NSEQ = 0;
    char *seq = NULL, *name = NULL;
    FASTAFILE *IN = OpenFasta(fastafile);

    if (IN == NULL) {
        fprintf(stderr, "Error opening input file %s\n", fastafile);
        return 1;
    }

    while (ReadFasta(IN, &seq, &name, &L) && Nseq < N) {
        sequences[Nseq] = safemalloc(sizeof(FASTASEQ));
        sequences[Nseq]->header = strdup(name);
        sequences[Nseq]->sequence = safemalloc((L+1)*sizeof(char));
        strncpy(sequences[Nseq]->sequence, seq, L);
        sequences[Nseq]->sequence[L] = '\0';
        sequences[Nseq]->length = L;

        safefree((void**)&seq);
        safefree((void**)&name);
        Nseq++;
    }

    CloseFasta(IN);
}

/**
 * Check if a string contains the character 'X'.
 *
 * @param motif The string to check.
 * @param size The size of the string.
 * @return 0 if no 'X' is found, -1 otherwise.
 */
int VerifyX(char * motif, int size) {
    for(int i = 0; i < size; i++) {
        if(motif[i] == 'X'){
            return -1;
        }
    }
    return 0;
}

/**
 * Create a copy of a string of a specified size.
 *
 * @param str The string to copy.
 * @param size The size of the string.
 * @return A pointer to the copied string.
 */
char * CopyString(char * str, int size) {
    char *newstr = safemalloc((size + 1) * sizeof(char));
    strncpy(newstr, str, size);
    newstr[size] = '\0';
    return newstr;
}

/**
 * Duplicate a string to an identical copy.
 *
 * @param str The string to copy.
 * @return A pointer to the copied string.
 */

char* DuplicateString(const char *str) {
/* Copy the Character Array given in argument in a new array with the same size (Character Array Duplication) */ 
	int n = strlen(str);
	char *dup = safemalloc((n+1)*sizeof(char));
	if(dup) {    strcpy(dup, str);    }
	dup[n]='\0';
	return dup;
}

/**
 * Compare two motifs of the same size.
 *
 * @param str1 The first motif to compare.
 * @param str2 The second motif to compare.
 * @return 0 if the motifs are equal or if str2 is '.', -1 if the lengths do not match, 1 otherwise.
 */
int CompareMotifs(char * str1, char * str2) {
    int S1 = strlen(str1); 
    int S2 = strlen(str2);

    if(S1 != S2) { 
        return -1;
    } else {
        for(int i = 0; i < S1; i++) {
            if(str1[i] != str2[i] && str2[i] != '.') {
                return 1;
            }
        }
        return 0;
    }
}

/**
 * Format a long integer to a string with leading zeros.
 *
 * @param N The width of the formatted string.
 * @param string The formatted string.
 * @param numb The number to format.
 */
void FormatLongInt(int N, char *string, long int numb) {
    
    if (N > MOTIFMAXLEN) {
        fprintf(stderr,"Motif size too big (>%d positions)\n",MOTIFMAXLEN); 
        exit(0);
    }
    char format[5];
    sprintf(format, "%%0%ldld", numb);
    sprintf(string, format, numb);
}

/**
 * Create degenerate versions of a motif.
 *
 * @param motif The motif to degenerate.
 * @param combinations The combinations to use for degeneration.
 * @param MOTIFLEN The size of the motif.
 * @param NC The number of combinations.
 * @return An array of degenerate motifs.
 */
char **DegenerateMotif(char *motif, char **combinations, int MOTIFLEN, int NC) {
    char ** degenerated = safemalloc(NC * sizeof(*degenerated));

    for(int comb = 0; comb < NC; comb++) {
        degenerated[comb] = safemalloc((MOTIFLEN + 1) * sizeof(**degenerated));
        for(int pos = 0; pos < MOTIFLEN; pos++) {
            degenerated[comb][pos] = (combinations[comb][pos] == '0') ? '.' : motif[pos];
        }
        degenerated[comb][MOTIFLEN] = '\0';
    }

    return degenerated;
}

/**
 * Convert a long integer to a binary string of a specified size.
 *
 * @param num The number to convert.
 * @param N The width of the binary string.
 * @return A pointer to the binary string.
 */
char * ToBinary(long int num, int N) {
    long int binary = 0, place = 0;
    char *integer = safemalloc((N + 1) * sizeof(char));
    integer[N] = '\0';  // Corrected this from N+1 to N.
    
    while(num != 0) {
        binary = binary + (num % 2 * pow(10, place));
        num = num / 2;
        place = place + 1;
    }

    FormatLongInt(N, integer, binary);
    return integer;
}

/**
 * Generate masks for degenerate motifs.
 *
 * @param size The size of the motifs.
 * @return A 2D array of masks.
 */
char ** MaskMotifPositions(int size) {
    int NC = 0;
    char **masks = NULL;
    
    int lastmask  = pow(2, size)          // last mask have 1 at all positions (i.e. defined position)
    int firstmask = (( lastmask / 2) + 3) // first mask must start with 1 and end with 1

    // Count the number of combinations (NC).
    for(int i = firstmask ; i < lastmask; i += 2) {
        NC++;
    }

    NC = (lastmask - firstmask) / 2; // Number of degenerate masks

    
    masks = safemalloc((NC - 1) * sizeof(*masks));
    
    // Generate the combinations.
    NC = 0;
    for(int i = firstmask; i < lastmask; i += 2) {
        masks[NC] = ToBinary((long int)i, size);
        NC++;
    }
    
    return masks;
}

/**
 * Scan a set of sequences for occurrences of a given motif of a specified size.
 *
 * @param k_Array A Judy array for k-mer occurrences.
 * @param m_Array A pointer to a Judy array for motif occurrences.
 * @param seqfile The sequence file name.
 * @param size The size of the motif.
 * @param NSEQ The number of sequences.
 * @param outfile The output file name.
 * @param sequences A pointer to the SEQUENCES structure.
 * @param AA The amino acid alphabet used.
 * @return 0 on successful completion.
 */
int ScanMotif(Pvoid_t k_Array, Pvoid_t *m_Array, char *seqfile, int size, int NSEQ, char *outfile, SEQUENCES *sequences, char *alphabet) {
    int s = 0, pos = 0, nletters = 0, firstpos = 0, lastpos = 0, NC = 0, mot = 0, i = 0;
    long int firstmask = 0, lastmask = 0;
    char *motif = NULL;
    char **masks = NULL;
    char **motifs = NULL;
    FILE* output = (FILE*)NULL;
    Pvoid_t tmpArray = (PWord_t)NULL;
    Pvoid_t occArray = (PWord_t)NULL;
    PWord_t occ, tmp, k, m;
    Word_t Rc_word;

	  nletters = 
    masks = MaskMotifPositions(size);
		
    //for(i = firstmask; i < lastmask; i += 2) { NC++; }

    for(s = 0; s < NSEQ; s++) {
        for(pos = 0; pos < sequences[s]->length - size - 1; pos++) {
            motif = CopyString(sequences[s]->sequence + pos, size);
            if( VerifyX(motif, size) ) {
                motifs = DegenerateMotif(motif, masks, size, NC);
                for(mot = 0; mot < NC; mot++) {
                    strcpy((char *)(Index), motifs[mot]);
                    Index[size] = '\0';
                    JSLG(k, k_Array, Index);
                    if(k != NULL) {
                        JSLG(tmp, tmpArray, Index);
                        JSLG(m, *m_Array, Index);
                        JSLG(occ, occArray, Index);
                        if(tmp != NULL) {
                            (*tmp)++;
                            (*occ)++;
                        } else {
                            JSLI(tmp, tmpArray, Index);
                            (*tmp)++;
                            if(occ != NULL && m != NULL) {
                                (*occ)++;
                                (*m)++;
                            } else {
                                JSLI(m, *m_Array, Index);
                                JSLI(occ, occArray, Index);
                                (*m)++;
                                (*occ)++;
                            }
                        }
                    }
                    safefree((void**)&motifs[mot]);
                }
                safefree((void**)&motifs);
            }
            safefree((void**)&motif);
        }
        JSLFA(Rc_word, tmpArray);
    }

    output = fopen(outfile, "a");
    if (!output) {
        fprintf(stderr, "Error! Unable to open the results file \"%s\".\n", outfile);
        exit(1);
    }


    for(firstpos=0;firstpos<nletters;firstpos++) {
		for(lastpos=0;lastpos<nletters;lastpos++) {
			Index[0] = '\0';
			JSLF(m, *(m_Array), Index);	/* get first string */
			while (m != NULL) {
				JSLG(occ,occArray,Index); 
				if(AA[firstpos] == Index[0] && AA[lastpos] == Index[size-1]){
					fprintf(output,"%s\t%d\t%d\n", Index,(int)(*occ),(int)(*m));
				}
				JSLN(m, *(m_Array), Index);   /* get next string */
			}
		}
	}
	
	safefree((void**)&combis);
	fclose(output);

	JSLFA(Rc_word, k_Array);		/* free array */
	JSLFA(Rc_word, occArray);		/* free array */
	return 0;
}

int ReadEnumeration(char *enumfile, char *delim, int s, Pvoid_t * success_Array, int onevsall) {

	FILE* enumfile = (FILE*)NULL;
	char line [MAXLINE]; char *pcourant=NULL;
	PWord_t success;     /* Judy array element. */

	ENUM = fopen(enumfile, "r");
	/* Check if the file could be open and read */
	if (!ENUM)	{	fprintf(stderr,"Error! Unable to open the delimited file of enumerated motifs \"%s\".\n",enumfile);	exit(1);	}

	while(fgets(line,MAXHEADER,ENUM)) {
		if( line[0] != '#' ){

			p = strtok(line, delim);
			if (p != NULL){
				strcpy((char *)(Index),p); Index[s]='\0';
			}

			p = strtok(NULL, delim);
			/* All occurrences per sequence */
			if(onevsall == 0){
				if (p != NULL){
					JSLI(success, *(success_Array),Index);   /* store string into array */
					(*success)=atoi(p);
				}
			}

			p = strtok(NULL, "\n");
			/* One occurrence per sequence */
			if(onevsall != 0){
				if (p != NULL){
					JSLI(success, *(success_Array),Index);   /* store string into array */
					(*success)=atoi(p);
				}
			}

		}
	}
	fclose(ENUM);
	return 0;
}

int main (int argc,char **argv){

	/* VARIABLES DECLARATION */
	int nseq; int size; int i; int onevsall;
	char * seqfile, *enumfile, *output;
	SEQUENCES * sequences;
	Pvoid_t  P_mArray, S_kArray; 

	/* INITIALIZE DEFAULT VALUE */
	NSEQ=0; size=0; i=0; onevsall=0;
	seqfile = NULL ; enumfile = NULL; output = NULL;
	sequences=NULL;
	P_mArray = (PWord_t) NULL;  /* Judy array. */
	S_kArray = (PWord_t) NULL;  /* Judy array. */

	/* READING INPUT */
	seqfile  =  strdup(argv[1]);
	enumfile =  strdup(argv[2]); 
	nseq     =  atoi(argv[3]); 
	size     =  atoi(argv[4]);
	output   =  strdup(argv[5]);
	if(argc == 7){
		onevsall = atoi(argv[6]);
	}
/*	printf("\n");*/
/*	fprintf(stderr,"INPUT :\n motif from set => %s \n motif from background => %s \n",enumfile,seqfile);*/

	sequences = safemalloc(NSEQ*sizeof(SEQUENCES));
	LoadFasta(seqfile,sequences,nseq);
	ReadEnumeration(enumfile,"\t",size,&S_kArray,onevsall);
	/* Scan each motif from the set in the Proteome and write a Motif enumeration file for Proteome */
	ScanMotif(S_kArray,&P_mArray,seqfile,size,nseq,output,sequences);

	for(i=0;i<nseq;i++){
		safefree((void**)&sequences[i]->header);
		safefree((void**)&sequences[i]->sequence);
		safefree((void**)&sequences[i]);
	}
	safefree((void**)&sequences);

	safefree((void**)&seqfile);
	safefree((void**)&enumfile); 
	safefree((void**)&output); 

/*	printf("\n");*/
	return EXIT_SUCCESS;
}
