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
/*                                   Motif Scan in Sequences                                  */
/*============================================================================================*/
/* Takes a motif enumeration file in input and scan those motifs in a set of sequences        */


/* Compile the following script */
/* gcc  ScanMotif.c safealloc.c fasta.c motif.c -o ScanMotif.exe -lm -lJudy -lpthread */ 

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
#include "fasta.h"

#define MIN(a,b)   (a < b ? a : b)
#define MAX(a,b)   (a > b ? a : b)
#define JUDYERROR_SAMPLE 1
#define LINELENGTH 50

uint8_t   Index[MAXLINE];   /* string to insert */

typedef struct {
  char* enumfile;
  char* fastafile;
  char* output;
  bool one;
  int verbose;
} OptionValues;

/* Main functions */
int ScanMotif(Pattern *, Pvoid_t, Pvoid_t *, SEQUENCES *, int, FILE *);  /*  Scan the proteome with a list of motifs, and count their occurences in the Proteome        */

void printUsage() {
  printf("Usage: MotifEnrichment.exe -e <enumfile> -f <fastafile> -o <output> [-g]\n");
  printf("Arguments (* = required):\n");
  printf("* -e, --enumfile  <enumfile>  Path to enumerated motifs from file\n");
  printf("* -f, --fastafile <fastafile> Path to FASTA formatted sequences file\n");
  printf("* -o, --output    <output>    Path to output file\n");
  printf("  -g, --one       <one>       Count one occurrence per sequence or all occurrences\n");
}

int checkOptionValues(const struct option* long_options, OptionValues* options_values) {
  if (options_values->enumfile == NULL) {
    fprintf(stderr, "Missing enumerated motifs!\n");
  }else{
    fprintf(stderr, "enumfile is set to  : %s\n", options_values->enumfile);
  }

  if (options_values->fastafile == NULL) {
    fprintf(stderr, "Missing fastafile!\n");
  }else{
    fprintf(stderr, "sequence is set to  : %s\n", options_values->fastafile);
  }

  if (options_values->output == NULL) {
    fprintf(stderr, "Missing output file!\n");
  }else{
    fprintf(stderr, "output is set to    : %s\n", options_values->output);
  }


  /* Check if required options are provided */
  if (options_values->enumfile == NULL || options_values->fastafile == NULL || options_values->output == NULL) {
    printUsage();
    exit(404);
  }

  return 0;
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
int ScanMotif(Pattern * pattern, Pvoid_t k_Array, Pvoid_t *m_Array, SEQUENCES *sequences, int nseq, FILE * OUT) {
  int s = 0, i=0, pos = 0, nletters = 0, p1=0, p2=0, first = 0, last = 0, opt = 0; int mlen=16;
  char *motif = NULL;
  char **motifs = NULL;
  
  Pvoid_t tmpArray = (PWord_t)NULL;
  Pvoid_t occArray = (PWord_t)NULL;
  PWord_t occ, tmp, k, m;
  Word_t Rc_word;

  Matrix masks = getMotifMask(pattern->len,pattern->cmin,pattern->cmax);
  
  for(s = 0; s < nseq; s++) {
    for(pos = 0; pos < sequences[s]->length - pattern->len - 1; pos++) {
      motif = CopyString(sequences[s]->sequence + pos, pattern->len);
      if( VerifyX(motif, pattern->len) == 0 ) {
        motifs = DegenerateMotif(masks.matrix, motif, pattern->len, masks.rows);
        for(i = 0; i <  masks.rows; i++) {
          strcpy((char*)(Index), motifs[i]);
          Index[mlen] = '\0';
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
          safefree((void**)&motifs[i]);
        }
        safefree((void**)&motifs);
      }
      safefree((void**)&motif);
    }
    JSLFA(Rc_word, tmpArray);
  }

  fprintf(stderr,"Scanned enumerated motifs in all sequences...\n");

  nletters = strlen(alphabet.NT);
  for(p1=0;p1<nletters;p1++) {
    first = alphabet.NT[p1];
    for(p2=0;p2<nletters;p2++) {
      last= alphabet.NT[p2];
      Index[0] = '\0';
      JSLF(m, *(m_Array), Index);  /* get first string */
      while (m != NULL) {
        JSLG(occ,occArray,Index);
        if( first == Index[0] && last == Index[pattern->len-1]){
          fprintf(OUT,"%s\t%d\t%d\n", Index,(int)(*occ),(int)(*m));
        }
        JSLN(m, *(m_Array), Index);   /* get next string */
      }
    }
  }
  
  freeMatrix(masks);
  JSLFA(Rc_word, k_Array);  /* free array */
  JSLFA(Rc_word, occArray);  /* free array */
  return 0;
}

int main (int argc,char **argv){

  /* VARIABLES DECLARATION AND INITIALIZATIONS*/
  int nseq = 0; int i = 0;  int verbose=0; int opt;
  char * fasta  = NULL;
  SEQUENCES * sequences = NULL;
  Pvoid_t  P_mArray = (PWord_t) NULL; 
  Pvoid_t  S_kArray = (PWord_t) NULL; 
  FILE * OUT = NULL;
  Pattern * motif = NULL;

  /* Variables to store the option values */
  OptionValues options = { NULL, NULL, NULL, false , 0};

  const char* short_options = "e:f:o:vg";
  /* Define the long options array */
  const struct option long_options[] = {
  { "enumfile", required_argument, NULL, 'e' },
  { "fastafile", required_argument, NULL, 'f' },
  { "output", required_argument, NULL, 'o' },
  { "one", no_argument, NULL, 'g'},
  { "verbose", no_argument, NULL, 'v'},
  { NULL, 0, NULL, 0 }
  };


  printf("-------------------INPUT-------------------\n");
  options.one = true;
  while ((opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
    switch (opt) {
    case 'e':
      options.enumfile = strdup(optarg);
      break;
    case 'f':
      options.fastafile = strdup(optarg);
      break;
    case 'o':
      options.output = strdup(optarg);
      break;
    case 'v':
      verbose=1;
      fprintf(stderr,"Verbose mode is enabled.\n");
      break;
    case 'g':
      options.one = true;
    case '?':
      printf("Unknown option: %c\n", optopt);
      exit(EXIT_FAILURE);
    default:
      printUsage();
      return 1;
    }
  }
  checkOptionValues(long_options, &options);
  printf("---------------------------------------------\n");

  // Process non-option arguments (if any)
  for (i = optind; i < argc; i++) {
   printf("Non-option argument: %s\n", argv[i]);
  }

  OUT = fopen(options.output, "w");
  if (OUT == NULL) {
    fprintf(stderr,"Error opening output file.\n");
    return -1;
  }
  
  nseq = countSequences(options.fastafile);
  sequences = safemalloc(nseq*sizeof(SEQUENCES));
  LoadFasta(options.fastafile,sequences,nseq);
  motif = ReadEnumeration(options.enumfile,"\t",&S_kArray,options.one);
  
  fprintf(OUT,"#Motifs        : %s\n",options.enumfile);
  fprintf(OUT,"#Input         : %s\n",options.fastafile);
  fprintf(OUT,"#Sequences     : %d\n",nseq);  
  fprintf(OUT,"#Pattern L/w/W : %d/%d/%d\n",motif->len,motif->xmin,motif->xmax);
  
  fprintf(stderr,"Scanned motif composition: %d X[%d,%d] not-X[%d,%d]\n",motif->len,motif->xmin,motif->xmax,motif->cmin,motif->cmax);

  /* Scan each motif from the set in the Proteome and write a Motif enumeration file for Proteome */
  ScanMotif(motif, S_kArray,&P_mArray,sequences,nseq,OUT);
  ClearFasta(sequences, nseq);

  safefree((void**)&options.enumfile); 
  safefree((void**)&options.fastafile); 
  safefree((void**)&options.output); 
  return EXIT_SUCCESS;
}
