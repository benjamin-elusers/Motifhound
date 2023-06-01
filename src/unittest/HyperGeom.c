/* +------------------------------------------------------------+ */
/* |  Date        : 10/05/2012                                  | */
/* |  Last Update : 28/05/2023                                  | */
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
/* gcc  MotifEnrichment.c safealloc.c fasta.c motif.c -o MotifEnrichment.exe -lm -lJudy -lpthread */ 

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "Judy.h"
#include "safealloc.h"
#include "motif.h"


#define MAXLINE 512 /* max string (line) length */
#define MAXHEADER 512 /* max string (line) length */
#define JUDYERROR_SAMPLE 1
#define LINELENGTH 50

uint8_t   Index[MAXLINE];   /* string to insert */


void printUsage() {
  printf("Usage: MotifEnrichment.exe -k <integer> -n <integer> -K <integer> -N <integer>\n");
  printf("Arguments (* = required):\n");
  printf("* -k, --set         <integer>  Occurrences in set\n");
  printf("* -n, --background  <integer>  Number of sequences in set\n");
  printf("* -K, --output      <integer>  Occurrences in background\n");
  printf("* -N, --nset        <integer>  Number of sequences in background\n");
}

int main (int argc,char **argv){

  int k = 0; int n = 0; int K=0; int N=0;
  int opt;
  long double pv = 0.0f;

  const char* short_options = "k:n:K:N:";
  /* Define the long options array */
  const struct option long_options[] = {
    { "k", required_argument, NULL, 'k' },
    { "n", required_argument, NULL, 'n' },
    { "K", required_argument, NULL, 'K' },
    { "N", required_argument, NULL, 'N' },
    { NULL, 0, NULL, 0 }
  };
  while ((opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
    switch (opt) {
    case 'k':
      k = atoi(optarg);
      break;
    case 'n':
      n = atoi(optarg);
      break;
    case 'K':
      K = atoi(optarg);
      break;
    case 'N':
      N = atoi(optarg);
      break;
    case '?':
      printf("Unknown option: %c\n", optopt);
      exit(EXIT_FAILURE);
    default:
      printUsage();
      return 1;
    }
  }
  fprintf(stderr,"k=%d\tx=%d\tm=%d\tn=%d\t\n",k,n,K,N);


  /* DECLARATION  AND INITIALIZATION*/
  pv = upperHG(k,n,K,N);
  fprintf(stderr,"k=%d\tx=%d\tm=%d\tn=%d\tpv=%.16Le (%.16Lf)\n",k,n,K,N,pv,pv);

  return EXIT_SUCCESS;
}

