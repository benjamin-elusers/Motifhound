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
/*                                  Motif Scan in Sequences                                   */
/*============================================================================================*/
/* Takes a motif enumeration file in input and scan those motifs in a set of sequences        */


/* Compile the following script */
/* gcc -Wall -ansi -pedantic -g ScanMotif.c -o ScanMotif.exe -lm -lJudy -lfasta -lsafealloc -lmotif -lpthread */ 

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <unistd.h> 
#include <getopt.h>
#include "Judy.h"
#include "safealloc.h"
#include "motif.h"
#include "fasta.h"

#define JUDYERROR_SAMPLE 1
#define MIN2(a,b)   (a < b ? a : b)

uint8_t   Index[MAXLINE];		 /* string to insert */

/* Main functions */
int    ScanMotif(Pvoid_t,Pvoid_t*,char*,int,int,char*,SEQUENCES*);	/*	Scan the proteome with a list of motifs, and count their occurences in the Proteome							*/
int    ReadEnumeration(char*,char*,int,Pvoid_t *,int);			/*	Read and store the enumeration file in a Judy Array															*/


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
int ScanMotif(Pvoid_t k_Array, Pvoid_t *m_Array, char *seqfile, int size, int NSEQ, char *outfile, SEQUENCES *sequences) {
    int s = 0, pos = 0, nletters = 0, firstpos = 0, lastpos = 0, NC = 0, mot = 0;
    char *motif = NULL;
    char **masks = NULL;
    char **motifs = NULL;
    FILE* output = (FILE*)NULL;
    Pvoid_t tmpArray = (PWord_t)NULL;
    Pvoid_t occArray = (PWord_t)NULL;
    PWord_t occ, tmp, k, m;
    Word_t Rc_word;

    masks = MaskMotifPositions(size);
		
    /* for(i = firstmask; i < lastmask; i += 2) { NC++; } */

    for(s = 0; s < NSEQ; s++) {
        for(pos = 0; pos < sequences[s]->length - size - 1; pos++) {
            motif = CopyString(sequences[s]->sequence + pos, size);
            if( VerifyX(motif, size) ) {
                motifs = DegenerateMotif(masks,motif, size, NC);
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
                if(alphabet[firstpos] == Index[0] && alphabet[lastpos] == Index[size-1]){
                    fprintf(output,"%s\t%d\t%d\n", Index,(int)(*occ),(int)(*m));
                }
                JSLN(m, *(m_Array), Index);   /* get next string */
            }
        }
    }
	
    safefree((void**)&masks);
    fclose(output);

    JSLFA(Rc_word, k_Array);		/* free array */
    JSLFA(Rc_word, occArray);		/* free array */
    return 0;
}

int ReadEnumeration(char *enumfile, char *delim, int s, Pvoid_t * success_Array, int onevsall) {

    FILE* ENUM = (FILE*)NULL;
    char line [MAXLINE]; char *p=NULL;
    PWord_t success;     /* Judy array element. */

    ENUM = fopen(enumfile, "r");
    /* Check if the file could be open and read */
    if (!ENUM)	{	fprintf(stderr,"Error! Unable to open the delimited file of enumerated motifs \"%s\".\n",enumfile);	exit(1);	}

    while(fgets(line,HEADERMAXLEN,ENUM)) {
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
	nseq=0; size=0; i=0; onevsall=0;
	seqfile = NULL ; enumfile = NULL; output = NULL;
	sequences=NULL;
	P_mArray = (PWord_t) NULL;  /* Judy array. */
	S_kArray = (PWord_t) NULL;  /* Judy array. */

	/* READING INPUT */
	seqfile  =  strdup(argv[1]);
	enumfile =  strdup(argv[2]); 
	size     =  atoi(argv[4]);
	output   =  strdup(argv[5]);
	if(argc == 7){
		onevsall = atoi(argv[6]);
	}

/*	printf("\n");*/
/*	fprintf(stderr,"INPUT :\n motif from set => %s \n motif from background => %s \n",enumfile,seqfile);*/

	sequences = safemalloc(nseq*sizeof(SEQUENCES));
	nseq = countSequences(seqfile);
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
