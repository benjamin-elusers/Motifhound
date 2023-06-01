/* fasta.h */
#ifndef FASTA_H
#define FASTA_H

#define FASTALINEWIDTH 80 /* max number of symbols in fasta formatted sequence */
#define HEADERMAXLEN 512 /* max header length */
#define SEQMAXLEN 40000 /* max sequence length */
#define MAXLINE 512 /* max string length */

typedef struct {
    const char AA[22+1]; 
    const char NT[6+1]; 
} ALPHABET;


typedef struct {
  FILE* fp;
  char  buffer[SEQMAXLEN];
} FASTAFILE;

typedef struct {
	char* header;   /* description/name */
	char* type;     /* NT or DNA or RNA or AA */
	char* sequence; /* symbols */
	int length;     /* length */
} SEQUENCE;

typedef SEQUENCE* SEQUENCES;

/* PROTOTYPES */

char* getAlphabet(const char *);
int   countSequences(const char*);               /* Count number of sequences in FASTA file    */
FASTAFILE* OpenFasta(const char*);               /* Open FASTA file for reading                */
int   ReadFasta(FASTAFILE*,char**,char**,int*);  /* Read sequences from FASTA file             */
int   WriteFasta(const char*,SEQUENCES*,int);    /* Write sequences to output FASTA file       */
int   LoadFasta(char*,SEQUENCES*,int);           /* Load FASTA sequences into SEQUENCES object */
void  CloseFasta(FASTAFILE*);                    /* Close FASTA file                           */
void  ClearFasta(SEQUENCES *, int);                   /* Clear sequences from memory                */

extern ALPHABET alphabet;

#endif