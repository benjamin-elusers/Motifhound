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
/* gcc  Enumeration.c safealloc.c fasta.c motif.c -o Enumeration.exe -lm -lJudy -lpthread */ 

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
  int mlen;
  int kmin;
  int xmin;
  int xmax;
  bool one;
} OptionValues;


#define MAXLINE 512 /* max string (line) length */
#define MAXHEADER 512 /* max string (line) length */
#define JUDYERROR_SAMPLE 1
#define LINELENGTH 50

uint8_t   Index[MAXLINE];   /* string to insert */

int   MotifEnrichment(double, Pvoid_t, Pvoid_t, char *, int, int);

void printUsage() {
  printf("Usage: MotifEnrichment.exe -s <set> -p <background> -o <output> -m <mlen> -k <kmin> [-x1 <xmin>] [-x2 <xmax>]\n");
  printf("Arguments (* = required):\n");
  printf("* -s, --set         <set>         Path to enumerated motifs from set file\n");
  printf("* -p, --background  <background>  Path to enumerated motifs from background file\n");
  printf("* -o, --output      <output>      Path to enumerated motifs from background file\n");
}

int checkOptionValues(const struct option* long_options, OptionValues* options_values) {
  for (int i = 0; long_options[i].name != NULL; i++) {
    if (options_values->set == NULL) {
      fprintf(stderr, "Missing enumerated motifs from set!\n");
      exit(404);
    }

    if (options_values->background == NULL) {
      fprintf(stderr, "Missing enumerated motifs from background!\n");
      exit(404);
    }

    if (options_values->output == NULL) {
      fprintf(stderr, "Missing output file!\n");
      exit(404);
    }
  }
    fprintf(stderr, "set is set to : %s\n", options_values->set);
    fprintf(stderr, "background is set to  : %s\n", options_values->background);
    fprintf(stderr, "output is set to   : %s\n", options_values->output);

    /* Check if required options are provided */
    if (options_values->set == NULL || options_values->background == NULL || options_values->output == NULL ||
    options_values->mlen == -1) {
    printUsage();
    exit(404);
    }

  return 0;
}


int MotifEnrichment(Pvoid_t k_Array, Pvoid_t m_Array,int nset,int nbg, FILE * OUT) {

  int w=0; 
  FILE* output = (FILE*) NULL;
  long double pv = 0.0f;
  PWord_t   k,m;


  /* Prints the motifs sorted lexicographically */
  Index[0] = '\0';
  JSLF(k, k_Array, Index);  /* get first string */
  while (k != NULL) {
  JSLG(m, m_Array, Index);
    if(m != NULL) {
      if((int)(*m) >= (int)(*k)){
        pv = upperHG((long int)(*k),(long int)(nset),(long int)(*m),(long int)(nbg));
        fprintf(OUT,"%s\t%d\t%d\t%d\t%d\t%.8Le\n", Index,(int)(*k),(int)(*m),(int)(nset),nbg,pv);
      }else{
        pv = upperHG((long int)(*k),(long int)(nset),(long int)(*m),(long int)(nbg));
        fprintf(OUT,"%s\t%d\t%d\t%d\t%d\t%.8Le | ( WARNING: This motif occurs more in the Set Sequences than in the Background Sequences )\n", Index,(int)(*k),(int)(*m),(int)(nset),nbg,pv);
      }
    }else{
      pv = upperHG((long int)(*k),(long int)(nset),(long int)(0+*k),(long int)(nbg));
      fprintf(OUT,"%s\t%d\t%d\t%d\t%d\t%.8Le | ( WARNING: This motif might occur less than 3 times in the Background Sequences )\n", Index,(int)(*k),(int)(*m),(int)(nset),nbg,pv);
    }
    JSLN(k, k_Array, Index);   /* get next string */
  }
  return 0;
}


int main (int argc,char **argv){

  /* DECLARATION  AND INITIALIZATION*/
  int nset = 0; int nbg = 0;
  Pvoid_t   P_mArray = (PWord_t)NULL;
  Pvoid_t   S_kArray = (PWord_t)NULL;
  Word_t    Rc_word;
  Motif * motif1 = NULL;
  Motif * motif2 = NULL;

  FILE* OUT = (FILE*)NULL;
  /* Variables to store the option values */
  OptionValues options = { NULL, NULL, NULL, -1, 3, 3, 3 };

  const char* short_options = "s:b:o:m:k:x:y:vg";
  /* Define the long options array */
  const struct option long_options[] = {
    { "set", required_argument, NULL, 's' },
    { "background", required_argument, NULL, 'b' },
    { "output", required_argument, NULL, 'a' },
    { "onevsall", no_argument, NULL, 'g'},
    { "verbose", no_argument, NULL, 'v'},
    { NULL, 0, NULL, 0 }
  };

  options.xmin = 3;
  options.kmin = 1;
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
  motif2 = ReadEnumeration(options.set,"\t",options.mlen,&S_kArray,options.one);

  OUT = fopen(options.output, "w");
  if (OUT == NULL) {
    fprintf(stderr,"Error opening output file.\n");
    exit(-1);
  }

/*  printf("INPUT :\n Sequences => %s \n Proteome => %s \n",Seqs_csvfile,Prot_csvfile);*/


  /* Calculate the hypergeometric cumulative distribution for each motif -> Pvalue */
  MotifEnrichment(options.mlen,S_kArray,P_mArray,nset,nbg,OUT);
  fclose(OUT);
  
  JSLFA(Rc_word,S_kArray);
  JSLFA(Rc_word,P_mArray);
  safefree(options.background);
  safefree(options.set); 
  return EXIT_SUCCESS;
}

