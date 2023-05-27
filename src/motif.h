/* motif.h */
#ifndef MOTIF_H
#define MOTIF_H
#define MOTIFMAXLEN 30
#define strdup(x) DuplicateString(x)

/* Functions for string manipulation */
char*  DuplicateString(const char*); /* Copy a string identically */
char*  CopyString(char*,int);        /* Copy the first n characters of a string */
void FormatLongInt(int, char *, unsigned long long int); /* format integer into a string format */
char*  ToBinary(unsigned long long int,int); /* Convert long integer into a binary format */


/* PROTOTYPES */
int    CompareMotifs(char *,char *); /* Compare two strings that may include wildcards (.) */
int    VerifyX(char*,int);           /* Check whether a motif contains the character 'X' */
char** DegenerateMotif(char**,char*,int,int); /* Degenerate motifs positions according to the masks */
char** MaskMotifPositions(int);               /* Create masks of degenerate motifs (i.e. with wildcards positions) */

#endif