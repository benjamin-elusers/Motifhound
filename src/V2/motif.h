/* motif.h */
#ifndef MOTIF_H
#define MOTIF_H
#include <stdbool.h>
#define MOTIFMAXLEN 30
#define MIN(a, b) (a < b ? a : b)
#define MAX(a, b) (a > b ? a : b)
#define strdup(x) DuplicateString(x)

typedef struct {
    char** matrix;
    int rows;
    int cols;
} Matrix;

typedef struct {
    int len;
    int cmin;
    int cmax;
    int xmin;
    int xmax;
} Pattern;


/* Functions for counting motifs */
int binomial(int, int);
long int longBinomial(long int, long int);
long double lowerHG(int, int, int, int);
long double upperHG(int, int, int, int);

/* Functions for string manipulation */
char* DuplicateString(const char*); /* Copy a string identically */
char* CopyString(char*, int);        /* Copy the first n characters of a string */
void FormatLongInt(int, char *, unsigned long long int); /* format integer into a string format */
char* ToBinary(unsigned long long int, int); /* Convert long integer into a binary format */


/* PROTOTYPES */
int CompareMotifs(char *, char *); /* Compare two strings that may include wildcards (.) */
int VerifyX(char*, int);           /* Check whether a motif contains the character 'X' */
int countX(char*,int);
char** DegenerateMotif(char**, char*, int, int); /* Degenerate motifs positions according to the masks */
char** MaskMotifPositions(int);               /* Create masks of degenerate motifs (i.e. with wildcards positions) */
Matrix createMatrix(int, int);
void freeMatrix(Matrix);
Matrix getMotifMask(int, int, int);
Pattern * ReadEnumeration(char *, char *, Pvoid_t *, bool);
char * ReadEnumerationParameters(char *enumfile);

#endif