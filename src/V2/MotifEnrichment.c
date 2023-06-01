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
#include "fasta.h"

typedef struct {
  char* set;
  char* background;
  char* output;
  int nset;
  int nbackground;
  bool one;
} OptionValues;


#define MAXLINE 512 /* max string (line) length */
#define MAXHEADER 512 /* max string (line) length */
#define JUDYERROR_SAMPLE 1
#define LINELENGTH 50

uint8_t   Index[MAXLINE];   /* string to insert */

int   MotifEnrichment(Pvoid_t, Pvoid_t, int, int, FILE *);

void printUsage() {
  printf("Usage: MotifEnrichment.exe -s <set> -p <background> -o <output> [-g]\n");
  printf("Arguments (* = required):\n");
  printf("* -s, --set         <filepath>         Path to enumerated motifs from set file\n");
  printf("* -b, --background  <filepath>  Path to enumerated motifs from background file\n");
  printf("* -o, --output      <filepath>  Path to enumerated motifs from background file\n");
  printf("* -n, --nset        <integer>   Number of sequences in the set file\n");
  printf("* -N, --nbackground <integer>   Number of sequences in the background file\n");
  printf("  -g, --one         <boolean>   Count one occurrence per sequence or all occurrences\n");
}

int checkOptionValues(const struct option* long_options, OptionValues* options_values) {
  if (options_values->set == NULL) {
    fprintf(stderr, "Missing enumerated motifs from set!\n");
  }else{
    fprintf(stderr, "set is set to         : %s\n", options_values->set);
  }

  if (options_values->background == NULL) {
    fprintf(stderr, "Missing enumerated motifs from background!\n");
  }else{
    fprintf(stderr, "background is set to  : %s\n", options_values->background);
  }

  if (options_values->output == NULL) {
    fprintf(stderr, "Missing output file!\n");
  }else{
    fprintf(stderr, "output is set to      : %s\n", options_values->output);
  }


  if (options_values->nset == -1) {
    fprintf(stderr, "Missing number of sequences in set file!\n");
  }else{
    fprintf(stderr, "nset is set to        : %d\n", options_values->nset);
  }

  if (options_values->nbackground == -1) {
    fprintf(stderr, "Missing number of sequences in background file!\n");
  }else{
    fprintf(stderr, "nbackground is set to : %d\n", options_values->nbackground);
  }

  
  /* Check if required options are provided */
  if (options_values->set == NULL || options_values->background == NULL || options_values->output == NULL || 
      options_values->nset == -1 || options_values->nbackground == -1 ) {
    printUsage();
    exit(404);
  }

  return 0;
}


int MotifEnrichment(Pvoid_t k_Array, Pvoid_t m_Array,int nset,int nbg, FILE * OUT) {

  int w=0; 
  int outlier=0;
  FILE* output = (FILE*) NULL;
  long double pv = 0.0f;
  PWord_t   k,m;


  /* Prints the motifs sorted lexicographically */
  Index[0] = '\0';
  JSLF(k, k_Array, Index);  /* get first string */
  while (k != NULL) {
  JSLG(m, m_Array, Index);
    if(m != NULL) {
      outlier = (int)(*k) >= (int)(*m) ? 1 : 0;
      if(!outlier){
        pv = upperHG((int)(*k),nset,(int)(*m),nbg);
        /* fprintf(stderr,"%s\t%d\t%d\t%d\t%d\t%.8Le\n", Index,(int)(*k),(int)(*m),(int)(nset),nbg,pv); */ 
        fprintf(OUT,"%s\t%d\t%d\t%d\t%d\t%.8Le\t%d\n", Index,(int)(*k),(int)(*m),(int)(nset),nbg,pv, outlier);
      }else{
        pv = upperHG((int)(*k),nset,(int)(*k)+(int)(*m),nset+nbg);
        fprintf(OUT,"%s\t%d\t%d\t%d\t%d\t%.8Le\t%d\n", Index,(int)(*k),(int)(*m),nset,nbg,pv, outlier);
      }
    }else{

      /* in R : phyper(k=nset,m=m,n=nbg-5,q=k-1,lower.tail = F) */
      pv = upperHG((int)(*k),nset,(int)(0+*k),nset+nbg);
      fprintf(OUT,"%s\t%d\t%d\t%d\t%d\t%.8Le\t%d\n", Index,(int)(*k),(int)(0),nset,nbg,pv,-1);
    }
    JSLN(k, k_Array, Index);   /* get next string */
  }
  return 0;
}


int main (int argc,char **argv){

  /* DECLARATION  AND INITIALIZATION*/
  int nset = 0; int nbg = 0; int len=0; int opt=0;
  int i = 0;
  int verbose;
  Pvoid_t   P_mArray = (PWord_t)NULL;
  Pvoid_t   S_kArray = (PWord_t)NULL;
  Word_t    Rc_word;
  Pattern * motif1 = NULL;
  Pattern * motif2 = NULL;

  FILE* OUT = (FILE*)NULL;
  /* Variables to store the option values */
  OptionValues options = { NULL, NULL, NULL, -1, -1, true };

  const char* short_options = "s:b:o:n:N:vg";
  /* Define the long options array */
  const struct option long_options[] = {
    { "set", required_argument, NULL, 's' },
    { "background", required_argument, NULL, 'b' },
    { "output", required_argument, NULL, 'a' },
    { "nset", required_argument, NULL, 'n' },
    { "nbackground", required_argument, NULL, 'N' },
    { "one", no_argument, NULL, 'g'},
    { "verbose", no_argument, NULL, 'v'},
    { NULL, 0, NULL, 0 }
  };

  verbose = 0;
  printf("-------------------INPUT-------------------\n");
  while ((opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
    switch (opt) {
    case 's':
      options.set = strdup(optarg);
      break;
    case 'b':
      options.background = strdup(optarg);
      break;
    case 'o':
      options.output = strdup(optarg);
      break;
    case 'n':
      options.nset = atoi(optarg);
      break;
    case 'N':
      options.nbackground = atoi(optarg);
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

  motif1 = ReadEnumeration(options.background,"\t",&P_mArray,options.one);
  motif2 = ReadEnumeration(options.set,"\t",&S_kArray,options.one);
  if(motif1->len != motif2->len || motif1->xmin != motif2->xmin || motif1->xmax != motif2->xmax){
    fprintf(stderr,"Error comparing non-identical patterns.\n");
    fprintf(stderr,"Pattern Background : L%d/%dX/%dX\n",motif1->len,motif1->xmin,motif1->xmax);
    fprintf(stderr,"Pattern Set        : L%d/%dX/%dX\n",motif2->len,motif2->xmin,motif2->xmax);
  }else{
    fprintf(stderr,"#Pattern    : L%d/%dX/%dX\n",motif1->len,motif1->xmin,motif1->xmax);  
  }
  
  OUT = fopen(options.output, "w");
  if (OUT == NULL) {
    fprintf(stderr,"Error opening output file.\n");
    exit(-1);
  }

  fprintf(OUT,"#Set        : %s\n",options.set);
  fprintf(OUT,"#Background : %s\n",options.background);  
  fprintf(OUT,"#Pattern    : L%d/%dX/%dX\n",motif1->len,motif1->xmin,motif1->xmax);  
  /* Calculate the hypergeometric cumulative distribution for each motif -> Pvalue */
  MotifEnrichment(S_kArray,P_mArray,options.nset,options.nbackground,OUT);
  fclose(OUT);
  
  JSLFA(Rc_word,S_kArray);
  JSLFA(Rc_word,P_mArray);
  safefree((void **)&options.background);
  safefree((void **)&options.set); 
  safefree((void **)&options.output); 

  return EXIT_SUCCESS;
}

