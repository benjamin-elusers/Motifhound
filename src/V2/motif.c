#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <stdbool.h>

#include "Judy.h"
#include "safealloc.h"
#include "motif.h"

#define JUDYERROR_SAMPLE 1
#define MOTIFMAXLEN 30
#define KEYLENGTH 50
#define MAXLINE 512

uint8_t Index[KEYLENGTH]; /* string to insert */

#define strdup(x) DuplicateString(x)

int binomial(int n, int k) {
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
    while (d > 1 && !(r % d)) {
      r /= d--;
    }
  }
  return r;
}

/*Binomial coefficient for large numbers*/
long int longBinomial(long int n, long int k) {
  long int r = 1;
  long int d = n - k;
  /* choose the smaller of k and n - k */
  if (d > k) {
    k = d;
    d = n - k;
  }

  while (n > k) {
    r *= n--;
    /* divide (n - k)! as soon as we can to delay overflows */
    while (d > 1 && !(r % d)) {
      r /= d--;
    }
  }
  return r;
}


/*Cumulative hypergeometric distribution*/
long double lowerHG(int k, int n, int m, int N) {
/* k: successes in the set, n: size of the set, m: successes in the population, N: size of the population */

  long double hg=1.0f; long double i; long double j;  long double bc=1.0f; /* long int  bc_tmp = 1.0f; */

  if(m <= N && k <= n && k <= m && n <= (N - m) && (n - k) <= (N - m)) {

    for (i=k; i<= 0; i-=1.0f) {
      bc = 1.0f;
      for (j=0.0f; j<n; j+=1.0f) {
        if (j < i)   bc *= (m - i + j + 1.0f) / (j + 1.0f);
        if (j < n-i) bc *= (N-m - n+i + j + 1.0f) / (j + 1.0f);
        if (j < n)   bc /= (N - n + j + 1.0f) / (j + 1.0f);
      }
      hg += bc;
      /* fprintf(stderr,"P(X=%d) = %.80Lf \t P(X>=%d) = %.80Lf \n",i,bc,i,hg); */
    }
/*    Return pvalue*/
    return (1.0f-hg);
  }else{
    return 0.0f;
  }
}

/*Cumulative hypergeometric distribution*/
long double upperHG(int k, int n, int m, int N) {
/* k: successes in the set, n: size of the set, m: successes in the population, N: size of the population */

  long double hg=0.0f; long double i; long double j; long double bc=1.0f;

  if(m <= N && k <= n && k <= m && n <= (N - m) && (n - k) <= (N - m)) {

    for (i=k; i<= MIN(n,m); i+=1.0f) {
      bc = 1.0f;
      for (j=0.0f; j<n; j+=1.0f) {
        if (j < i)   bc *= (m - i + j + 1.0f) / (j + 1.0f);
        if (j < n-i) bc *= (N-m - n+i + j + 1.0f) / (j + 1.0f);
        if (j < n)   bc /= (N - n + j + 1.0f) / (j + 1.0f);
      }
      hg += bc;
      /* fprintf(stderr,"P(X=%ld) = %.80Lf \t P(X>=%ld) = %.80Lf \n",i,bc,i,hg); */
    }
/*    Return pvalue*/
    return (hg);
  }else{
    return 0.0f;
  }
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
  if(dup) {  strcpy(dup, str);  }
  dup[n]='\0';
  return dup;
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
 * Format a long integer to a string with leading zeros.
 *
 * @param N The width of the formatted string.
 * @param string The formatted string.
 * @param numb The number to format.
 */
void FormatLongInt(int N, char *string, unsigned long long int numb) {
  
  char format[7];

  if (N > sizeof(long long int)) {
    fprintf(stderr,"Motif size too big (>%d positions)\n",MOTIFMAXLEN); 
    N = MOTIFMAXLEN;
  }
  sprintf(format, "%%0%dlld", N);
  sprintf(string, format, numb);
}


/**
 * Convert a long integer to a binary string of a specified size.
 *
 * @param num The number to convert.
 * @param N The width of the binary string.
 * @return A pointer to the binary string.
 */
char * ToBinary(unsigned long long int num, int N) {
  /*int place = 0;*/
  char *binary = safemalloc((N + 1) * sizeof(char));
  
  for (int i = N-1; i >= 0; i--) {
    binary[i] = (num & 1) + '0';
    num >>= 1;
  }
  binary[N] = '\0';


/*  while(num != 0) {
    binary = binary + (num % 2 * pow(10, place));
    num = num / 2;
    place = place + 1;
  }

  FormatLongInt(N, binary, num);
*/
  /* fprintf(stderr,"%s",binary); */
  return(binary);
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
  int i = 0;

  if(S1 != S2) { 
    return -1;
  } else {
    for(i=0; i < S1; i++) {
      if(str1[i] != str2[i] && str2[i] != '.') {
        return 1;
      }
    }
    return 0;
  }
}

/**
 * Check if a string contains the character 'X'.
 *
 * @param motif The string to check.
 * @param size The size of the string.
 * @return 0 if no 'X' is found, -1 otherwise.
 */
int VerifyX(char * motif, int size) {
  int i = 0;

  for(i=0; i < size; i++) {
    if(motif[i] == 'X'){
      return -1;
    }
  }
  return 0;
}



int countX(char * motif, int len){
  int p=0; int nX=0;
  
  for(p=1;p<len-1;p++){
    if(motif[p] == '.' || motif[p] == 'X' ){
     nX++; 
   }
  }
  
  return (nX);
}

/**
 * Create degenerate versions of a motif.
 *
 * @param combinations The combinations to use for degeneration.
 * @param motif The motif to degenerate.
 * @param len The size of the motif.
 * @param NC The number of combinations.
 * @return An array of degenerate motifs.
 */
char **DegenerateMotif(char **combinations, char *motif, int len, int NC) {
  int comb = 0, pos = 0;
  char ** degenerated = safemalloc(NC * sizeof(*degenerated));

  for(comb=0; comb < NC; comb++) {
    degenerated[comb] = safemalloc((len + 1) * sizeof(**degenerated));
    for(pos=0; pos < len; pos++) {
      degenerated[comb][pos] = (combinations[comb][pos] == '0') ? '.' : motif[pos];
    }
    degenerated[comb][len] = '\0';
  }

  return degenerated;
}


/**
 * Generate masks for degenerate motifs.
 *
 * @param size The size of the motifs.
 * @return A 2D array of masks.
 */
char ** MaskMotifPositions(int size) {
  int NC = 0, i = 0;
  int lastmask  = pow(2, size);      /* last mask have 1 at all positions (i.e. defined position) */
  int firstmask = (( lastmask / 2) + 3); /* first mask must start with 1 and end with 1 */
  char **masks = NULL;

  /* Count the number of combinations (NC). */
  for(i=firstmask; i < lastmask; i += 2) {
    NC++;
  }

  NC = (lastmask - firstmask) / 2; /* Number of degenerate masks */
  masks = safemalloc((NC - 1) * sizeof(*masks));
  
  /* Generate the combinations. */
  NC = 0;
  for(i=firstmask; i < lastmask; i += 2) {
    masks[NC] = ToBinary((unsigned long long int) i, size);
    NC++;
  }
  
  return masks;
}


Matrix createMatrix(int N1, int N2) {
  int i = 0;
  Matrix mat;
  
  mat.rows = N1;
  mat.cols = N2;
  // Dynamically allocate memory for the matrix
  mat.matrix = safemalloc(N1 * sizeof(char*));
  for (i=0; i < N1; i++) {
    mat.matrix[i] = safemalloc(N2 * sizeof(char));
  }
  return mat;
}

void freeMatrix(Matrix matrix) {
  for (int i = 0; i < matrix.rows; i++) {
    safefree((void**)&matrix.matrix[i]);
  }
  safefree((void**)&matrix.matrix);
}

Matrix getMotifMask(int mlen, int xmin, int xmax) {
  
  int i=0; 
  int j=0;
  int k=0;
  int x=0;

  int cmin = 0;
  int cmax = 0;
  int nc = 0;
  int nx = 0;
  int ncx = 0;
  int verbose = 0; int debug=0;

  cmin = pow(2, mlen - 1) + pow(2, 1) + pow(2, 0);
  cmax = cmin;
  for (i = 2; i < mlen - 1; i++) {
    cmax += pow(2, i);
  }
  nc = ((cmax - cmin) / 2) + 1;
  if(verbose){ 
    fprintf(stderr, "# Combinations %d \n", cmax - cmin);
    fprintf(stderr, "# Masks (exclude those with 0 in Cter) %d \n", nc); 
  }

  /* Compute number of masks with increasind number of fixed positions (i.e. decreasing wildcard positions) */ 
  for (i = xmin; i <= xmax && i <= mlen - 2; i++) {
    ncx += binomial(mlen - 2, i);
  }
  if(verbose){ fprintf(stderr, "# expected masks %d \n", ncx); }

  Matrix combinations = createMatrix(nc,mlen);
  Matrix masks = createMatrix(ncx,mlen);
  i = 0;
  for (j = cmin; j <= cmax; j += 2) {
    nx = 0;
    combinations.matrix[i] = ToBinary((unsigned long long int)(j), mlen);
    for (k = 0; k < mlen; k++) {
      nx += combinations.matrix[i][k] - '0';
    }
    if (mlen-nx >= xmin && mlen-nx <= xmax) {
      masks.matrix[x] = strdup(combinations.matrix[i]);
      x++;
    }
    i++;
  }

  if(verbose){ 
    fprintf(stderr, "# motif combinations %d\n", nc);
    fprintf(stderr, "from %s to %s\n", ToBinary(cmin, mlen), ToBinary(cmax, mlen));
    fprintf(stderr, "# motif combinations between %d and %d wildcard positions %d\n",xmin, xmax, ncx);
  }

  if(debug){
  	for(i=0; i<masks.rows; i++){
      fprintf(stderr, "VALID MASK %d. %s\n", i + 1, masks.matrix[i]);
    }
  }

  fprintf(stderr, "Number of %d-masks: %d\n", masks.cols, masks.rows); 
  fprintf(stderr, "Pattern length    : %d\n",mlen); 
  fprintf(stderr, "Wildcards         : [%d,%d]\n",xmin, xmax);
  fprintf(stderr, "Fixed             : [%d,%d]\n",mlen-xmax, mlen-xmin);
  freeMatrix(combinations);

  return(masks);
}

Pattern * ReadEnumeration(char *enumfile, char *delim, Pvoid_t * success_Array, bool one) {
  int len=0; int x1=999; int x2=0; int nx=0;
  FILE* ENUM = (FILE*)NULL;
  char line [MAXLINE]; 
  char *p=NULL;
  PWord_t success;
  Pattern * scanned = safemalloc(sizeof(Pattern));

  ENUM = fopen(enumfile, "r");
  /* Check if the file could be open and read */
  if (!ENUM)  {  fprintf(stderr,"Error! Unable to open the delimited file of enumerated motifs \"%s\".\n",enumfile);  exit(1);  }

  while(fgets(line,MAXLINE,ENUM)) {
    if( line[0] != '#' ){
      p = strtok(line, delim);
      if (p != NULL){
        len=strlen(p);
        strcpy((char *)(Index),p); 
        Index[len]='\0';
        nx = countX(Index,len);
        if(x1 > nx){ x1 = nx; }
        if(x2 < nx){ x2 = nx; }
      }

      p = strtok(NULL, delim);
      /* All occurrences per sequence */
      if(one == false){
        if (p != NULL){
          JSLI(success, *(success_Array),Index);   /* store string into array */
          (*success)=atoi(p);
        }
      }

      p = strtok(NULL, "\n");
      /* One occurrence per sequence */
      if(one == true){
        if (p != NULL){
          JSLI(success, *(success_Array),Index);   /* store string into array */
          (*success)=atoi(p);
        }
      }
    }
  }

  scanned->len  = len;
  scanned->xmin = x1;
  scanned->xmax = x2;
  scanned->cmin = len-x2;
  scanned->cmax = len-x1;

  fclose(ENUM);
  return (scanned);
}
