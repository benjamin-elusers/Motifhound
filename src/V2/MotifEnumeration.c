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
/*                             Motif Enumeration in Sequences                                 */
/*============================================================================================*/
/* Takes a motif enumeration file in input and scan those motifs in a set of sequences        */


/* Compile the following script */
/* gcc  MotifEnumeration.c safealloc.c fasta.c motif.c -o MotifEnumeration.exe -lm -lJudy -lpthread */ 

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "Judy.h"
#include "safealloc.h"
#include "motif.h"
#include "fasta.h"

#define JUDYERROR_SAMPLE 1
#define KEYLENGTH 50

uint8_t Index[KEYLENGTH]; /* string to insert */

typedef struct {
  char* fastafile;
  char* output;
  char* seqtype;
  int mlen;
  int kmin;
  int xmin;
  int xmax;
} OptionValues;

void printUsage();
int checkOptionValues(const struct option*, OptionValues*);
int sumArray(const int*, int);
int binomial(int, int);
int * Enumeration(int, Matrix, int, char, char, SEQUENCES*, int, FILE*);

void printUsage() {
  printf("Usage: Enumeration.exe -f <fastafile> -o <output> -m <mlen> -k <kmin> [-x1 <xmin>] [-x2 <xmax>] [-a <seqtype>]\n");
  printf("Arguments (* = required):\n");
  printf("* -f, --fasta   <fastafile>   Path to the input FASTA file\n");
  printf("* -o, --output  <outputfile>  Path to the output file\n");
  printf("* -a, --seqtype   <seqtype>   Sequence type (AA or NT)\n");
  printf("* -m, --mlen    <mlen>    Motif length\n");
  printf("  -k, --kmin    <kmin>    Minimal number of occurrences\n");
  printf("  -x, --xmin    <xmin>    Minimum number of wildcard positions\n");
  printf("  -y, --xmax    <xmax>    Maximum number of wildcard positions\n");
}

int checkOptionValues(const struct option* long_options, OptionValues* options_values) {
  if (options_values->fastafile == NULL) {
    fprintf(stderr, "Missing fastafile!\n");
  }else{
    fprintf(stderr, "fastafile is set to : %s\n", options_values->fastafile);
  }

  if (options_values->output == NULL) {
    fprintf(stderr, "Missing output file!\n");
  }else{
    fprintf(stderr, "output is set to    : %s\n", options_values->output);
  }

  if (options_values->seqtype == NULL) {
    fprintf(stderr, "Missing type of fasta sequences!\n");
  }else{
    fprintf(stderr, "seqtype is set to   : %s\n", options_values->seqtype);
  }

  if (options_values->mlen == -1) {
    fprintf(stderr, "Missing motif length!\n");
  }else{
    fprintf(stderr, "mlen is set to      : %d\n", options_values->mlen);
  }

  if (options_values->kmin == -1) {
    fprintf(stderr, "Missing minimum of motif occurrences (default to 1)!\n");
  }else{
    fprintf(stderr, "kmin is set to      : %d\n", options_values->kmin);
  }

  if (options_values->xmin < 0) {
    options_values->xmin = 0;
    fprintf(stderr, "Missing minimum wildcard positions (default to 0)!\n");
  }else{
    fprintf(stderr, "xmin is set to      : %d\n", options_values->xmin);
  }


  if (options_values->xmax < 0) {
    fprintf(stderr, "Missing maximum wildcard positions (default to motif size - 2)!\n");
    options_values->xmax = options_values->mlen -2;
  }else{
    fprintf(stderr, "xmax is set to      : %d\n", options_values->xmax);
  }


  /* Check if required options are provided */
  if (options_values->fastafile == NULL || options_values->seqtype == NULL || options_values->output == NULL ||
    options_values->mlen == -1) {
    printUsage();
    exit(404);
  }

  return 0;
}


int sumArray(const int* array, int length) {
  int sum = 0;
  for (int i = 0; i < length; i++) {
    sum += array[i];
  }
  return sum;
}

int * Enumeration(int len, Matrix masks, int K_OCC, char firstletter, char lastletter,
         SEQUENCES* sequences, int nseq, FILE* OUT) {
  int k = 0;
  int pos = 0;
  int s = 0;
  int * nmotifs = NULL;
  clock_t start, end;
  double cpu_time_used;

  char* motif = NULL;
  char** all_motifs = NULL;

  Pvoid_t occArray = (PWord_t)NULL;
  Pvoid_t seqArray = (PWord_t)NULL;
  Pvoid_t tmpArray = (PWord_t)NULL;

  PWord_t occ;
  PWord_t seq;
  PWord_t tmp;

  Word_t Rc_word;
  start = clock();  
  nmotifs = safecalloc(0,2 * sizeof(int));

  for (s = 0; s < nseq; s++) {
    for (pos = 0; pos < sequences[s]->length - len - 1; pos++) {
      if (sequences[s]->sequence[pos] == firstletter && sequences[s]->sequence[pos + len - 1] == lastletter) {
        motif = CopyString(sequences[s]->sequence + pos, len);
        if (VerifyX(motif, len) == 0) {
          all_motifs = DegenerateMotif(masks.matrix, motif, len, masks.rows);
          for (k = 0; k < masks.rows; k++) {
            strcpy((char*)(Index), all_motifs[k]);
            Index[len] = '\0';
            JSLG(tmp, tmpArray, Index);
            JSLG(occ, occArray, Index);
            JSLG(seq, seqArray, Index);
            if (tmp != NULL) {
              (*tmp)++;
              (*occ)++;
            } else {
              if (occ != NULL && seq != NULL) {
                (*occ)++;
                (*seq)++;
              } else {
                JSLI(seq, seqArray, Index); /* store string into array */
                (*seq)=1;
                JSLI(occ, occArray, Index); /* store string into array */
                (*occ)=1;
              }
              JSLI(tmp, tmpArray, Index); /* store string into array */
              (*tmp)=1;
            }
            safefree((void**)&all_motifs[k]); /* free array */
          }
          safefree((void**)&motif);
          safefree((void**)&all_motifs); /* free array */
        }
      }
    }
    JSLFA(Rc_word, tmpArray);
  }

  /* Prints the motifs sorted lexicographically */
  Index[0] = '\0';

  JSLF(seq, seqArray, Index); /* get first string */
  while (seq != NULL) {
    JSLG(occ, occArray, Index); /* get first string */
    if (occ != NULL) {
      nmotifs[0]++;
      if ((int)(*seq) >= K_OCC) {
        nmotifs[1]++;
        fprintf(OUT, "%s\t%d\t%d\n", Index, (int)(*occ), (int)(*seq));
      }
    }
    JSLN(seq, seqArray, Index); /* get next string */
  }
  JSLFA(Rc_word, seqArray); /* free array */
  JSLFA(Rc_word, occArray); /* free array */
  
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  if( nmotifs[0] > 0){
    fprintf(stderr, "%c-%c\t%12d\t%8d\t%3.1f\n",firstletter,lastletter,nmotifs[0],nmotifs[1],cpu_time_used);
  }

  return(nmotifs);
}


int main(int argc, char** argv) {
  /* DECLARATIONS & INITIALISATIONS */
  int i = 0;
  int j = 0;
  int first = 0;
  int last = 0;
  int opt;
  int nseq = 0;
  int npos = 0;
  int nletters = 0;
  int nmasks = 0;
  int verbose = 0;
  int rows = 0;
  int cols = 0;


  clock_t t0, t1;
  double total_time_used;
  int *nmotifs;
  int *NMOTIFS;

  char* letters = NULL;
  Matrix masks;
  SEQUENCES* sequences = NULL;
  FILE* OUT = (FILE*)NULL;

  /* Variables to store the option values */
  OptionValues options = { NULL, NULL, NULL, -1, 0, 0, 0};

  const char* short_options = "f:o:m:k:x:y:a:v";
  /* Define the long options array */
  const struct option long_options[] = {
    { "fasta", required_argument, NULL, 'f' },
    { "output", required_argument, NULL, 'o' },
    { "seqtype", required_argument, NULL, 'a' },
    { "mlen", required_argument, NULL, 'm' },
    { "kmin", optional_argument, NULL, 'k' },
    { "xmin", optional_argument, NULL, 'x' },
    { "xmax", optional_argument, NULL, 'y' },
    { "verbose", no_argument, NULL, 'v'},
    { NULL, 0, NULL, 0 }
  };

  verbose = 0;
  printf("-------------------INPUT-------------------\n");
  while ((opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
    switch (opt) {
    case 'f':
      options.fastafile = strdup(optarg);
      break;
    case 'o':
      options.output = strdup(optarg);
      break;
    case 'a':
      options.seqtype = strdup(optarg);
      break;
    case 'm':
      options.mlen = atoi(optarg);
      options.xmax = MIN(options.xmax,options.mlen - 2);
      break;
    case 'k':
      if (optarg != NULL) {
        options.kmin = atoi(optarg);
      }
      break;
    case 'v':
      verbose=1;
      fprintf(stderr,"Verbose mode is enabled.\n");
      break;
    case 'x':
      if (optarg != NULL) {
        options.xmin = MAX(0,atoi(optarg));
      }
      break;
    case 'y':
      if (optarg != NULL) {
        options.xmax = MIN(options.mlen-2,atoi(optarg));
      }
      break;
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

  letters = getAlphabet(options.seqtype);
  nletters = strlen(letters);
  if(verbose){ 
    fprintf(stderr, "Alphabet has %d letters : [%s]\n", nletters, letters);
  }

  nseq = countSequences(options.fastafile);
  sequences = safemalloc(nseq * sizeof(SEQUENCES));
  npos = LoadFasta(options.fastafile, sequences, nseq);

  fprintf(stderr, "Number of %d-mers : %d\n", options.mlen,npos);
  masks = getMotifMask(options.mlen, options.xmin, options.xmax);
  fprintf(stderr, "Motif space of length %d: %d (%d x %d)\n", options.mlen,npos*masks.rows,npos,masks.rows);


  printf("---------------------------------------------\n");
  t0 = clock();  
  printf("%3s\t%12s\t%8s(k>%d)\t%4s\n","S-E","# MOTIFS","# MOTIFS",options.kmin-1,"t(sec)");
  printf("=============================================\n");

  fprintf(OUT,"#Input     : %s\n",options.fastafile);  
  fprintf(OUT,"#Alphabet  : %s\n",options.seqtype);  
  fprintf(OUT,"#Sequences : %d\n",nseq);
  fprintf(OUT,"#Pattern   : L%d/%dX/%dX\n",options.mlen, options.xmin, options.xmax);  
  fprintf(OUT,"#Occurence : %d\n",options.kmin);  

  NMOTIFS=safecalloc(0,2);
  for (first = 0; first < nletters; first++) {
    for (last = 0; last < nletters; last++) {
      nmotifs = Enumeration(options.mlen, masks, options.kmin, letters[first], letters[last], sequences, nseq, OUT);
      NMOTIFS[0] += nmotifs[0];
      NMOTIFS[1] += nmotifs[1];
    }
  }

  printf("---------------------------------------------\n");
  t1 = clock();
  total_time_used = ((double) (t1 - t0)) / CLOCKS_PER_SEC;
  fprintf(stderr, "Total \t%12d\t%8d\t%3.1f\n",NMOTIFS[0],NMOTIFS[1],total_time_used);
  fprintf(OUT, "#MOTIFS (k>%d) %12d (%8d)\n",options.kmin,NMOTIFS[0],NMOTIFS[1]);
  fprintf(OUT, "#TIME %3.1f\n",total_time_used);
  fclose(OUT);

  ClearFasta(sequences,nseq);
  safefree((void**)&masks);
  safefree((void**)&options.seqtype);
  safefree((void**)&options.fastafile);
  safefree((void**)&options.output);

  return EXIT_SUCCESS;
}
