#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

#include "safealloc.h"
#include "motif.h"
#define MOTIFMAXLEN 30
#define strdup(x) DuplicateString(x)

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
	if(dup) {	strcpy(dup, str);	}
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
	binary[N+1] = '\0';


/*	while(num != 0) {
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
	int lastmask  = pow(2, size);		  /* last mask have 1 at all positions (i.e. defined position) */
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
