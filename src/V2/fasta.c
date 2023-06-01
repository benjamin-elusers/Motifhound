#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "safealloc.h"
#include "fasta.h"

#define FASTALINEWIDTH 80 /* max number of symbols in fasta formatted sequence */
#define HEADERMAXLEN 512 /* max header length */
#define SEQMAXLEN 40000 /* max sequence length */
#define MAXLINE 512 /* max string length */


ALPHABET alphabet = {
  .AA = "ACDEFGHIKLMNPQRSTUVWXY",
  .NT = "ACGNTU"
};

int countSequences(const char* filename) {
  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    printf("Error opening file: %s\n", filename);
    return -1;
  }

  int count = 0;
  char line[MAXLINE];
  while (fgets(line, sizeof(line), file) != NULL) {
    if (line[0] == '>') {
      count++;
    }
  }

  fclose(file);
  return count;
}

char * getAlphabet(const char * type){

  if (strcmp(type, "AA") == 0) {
    return("ACDEFGHIKLMNPQRSTUVWXY");
  } else if (strcmp(type, "NT") == 0) {
    return("ACGNTU");
  }else{
    fprintf(stderr,"Error unknown alphabet.\n");
    return NULL;
  }
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
  /* Allocate memory for the FASTAFILE structure */
  FASTAFILE* ffp = safemalloc(sizeof(FASTAFILE));

  /* Try to open the file */
  ffp->fp = fopen(seqfile, "r");

  /* If the file could not be opened, print an error message
     free the memory allocated for the FASTAFILE structure
     return NULL */

  if (ffp->fp == NULL){
    fprintf(stderr, "Error! Unable to open the FASTA formatted Sequences file \"%s\".\n", seqfile);
    safefree((void**)&ffp);
    return NULL;
  }

  /* Try to read the first line of the file
     If the line could not be read, 
     free the memory allocated for the FASTAFILE structure
     return NULL */
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
  /* Allocate string in blocks of SEQMAXLEN residues */
  char* sequence = safemalloc(sizeof(char) * SEQMAXLEN); 
  char* s = NULL; 
  char* description = NULL;
  int nalloc = SEQMAXLEN;
  int n = 0; 

  /* Exit if line buffer does not contain '>' */
  if (INPUT->buffer[0] != '>') {
    safefree((void**)&sequence);
    return 0; 
  }

  /* Store the description line */
  s = strtok(INPUT->buffer + 1, "\n");
  description = strdup(s);  

  /* Read the sequence */
  while(fgets(INPUT->buffer, FASTALINEWIDTH, INPUT->fp)){
    /* If the next description line was reached, stop reading */
    if (INPUT->buffer[0] == '>') break;
    
    /* Remove the newline character from the end of the buffer */
    INPUT->buffer[strlen(INPUT->buffer) - 1] = '\0';

    /* Store the sequence in an array of char (sequence) */
    for(s = INPUT->buffer; *s != '\0'; s++){ 
      sequence[n] = *s;
      n++;
      
      /* If there is no room left to store the sequence in the array of char, realloc */
      if (nalloc == n){
        nalloc += SEQMAXLEN;
        sequence = saferealloc(sequence, sizeof(char) * nalloc);
      }
    }
  }

  /* Set return values */
  *DESC = description;
  *SEQ = sequence;
  *L = n;

  return 1;
}

/**
 * Write sequences to a FASTA formatted file.
 *
 * This function writes sequences to a FASTA formatted file. Each sequence is written
 * with a line width of 80 (FASTALINEWIDTH) characters.
 *
 * @param fastafile The path of the FASTA file to write to.
 * @param sequences A pointer to an array of SEQUENCES structures.
 * @param N The number of sequences to write.
 * @return 0 on success, non-zero on failure.
 */
int WriteFasta(const char * fastafile, SEQUENCES *sequences, int N) {
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
int LoadFasta(char * fastafile, SEQUENCES * sequences, int N) {
  
  int L = 0, n = 0, np=0, nseq=0;
  char *seq = NULL, *name = NULL;
  FASTAFILE *IN = OpenFasta(fastafile);

  if (IN == NULL) {
    fprintf(stderr, "Error opening input file %s\n", fastafile);
    return 1;
  }

  nseq = countSequences(fastafile);
  fprintf(stderr, "Number of sequences found: %d\n", nseq);

  while (ReadFasta(IN, &seq, &name, &L) && n < N) {
    sequences[n] = safemalloc(sizeof(SEQUENCE));
    sequences[n]->header = strdup(name);
    sequences[n]->sequence = strdup(seq);
    sequences[n]->length = L;
    safefree((void**)&seq);
    safefree((void**)&name);
    np += L;
    n++;
  }

  CloseFasta(IN);
  return np;
}

void ClearFasta(SEQUENCES * sequences, int nseq){
  int i=0;
  for (i = 0; i < nseq; i++) {
  safefree((void**)&sequences[i]->header);
  safefree((void**)&sequences[i]->sequence);
  safefree((void**)&sequences[i]);
  }
  safefree((void**)&sequences);
}
