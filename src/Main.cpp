#include "Classe.hpp"
//---------------------------------------------------------------------------
// main routine

int main(int argc, char **argv){

	if (argc <= 1 ){
		printf("Usage : -pro proteome -pos positives -mot motifs\n\n");
		printf("-pro [input proteome (fasta file)] (required)\n");
		printf("-pos [input positives (fasta file)] (required)\n");
		printf("-mot [output motifs file] (required)\n");
		printf("-frq [minimum number of positives including the motif (defaule 3)] (optional)\n");
		printf("-min [minimum length of a motif (default 3)] (optional)\n");
		printf("-max [maximum length of a motif (default 8)] (optional)\n");
		printf("-nwc [minimum number of non wildcards in the motif (default 3)] (optional)\n");
		printf("-thr [number of threads (default 5)] (optional)\n");
		return(1);
	}else{
		char pro[MAXLEN];
		char pos[MAXLEN];
		char mot[MAXLEN];

		for (int i=0; i<argc; i++){
			if (strcmp(argv[i], "-pro") == 0){
				i++;
				// Get proteome file
				strcpy(pro, argv[i]);
			}else if(strcmp(argv[i], "-pos") == 0){
				i++;
				// Get positives file
				strcpy(pos, argv[i]);
			}else if (strcmp(argv[i], "-mot") == 0){
				i++;
				// Get motifs file
				strcpy(mot, argv[i]);
			}else if (strcmp(argv[i], "-frq") == 0){
				i++;
				// Get minimum number of positive including a motif
				FRQMOT = atoi(argv[i]);
			}else if (strcmp(argv[i], "-min") == 0){
				i++;
				// Get minimum length of a motif
				MINMOT = atoi(argv[i]);
			}else if (strcmp(argv[i], "-max") == 0){
				i++;
				// Get maximum length of a motif
				MAXMOT = atoi(argv[i]);
			}else if (strcmp(argv[i], "-nwc") == 0){
				i++;
				// Get minimum number of non-wild cards in a motif
				NWCMOT = atoi(argv[i]);
			}else if (strcmp(argv[i], "-thr") == 0){
				i++;
				// Get number of simultaneous threads
				THREAD = atoi(argv[i]);
			}
		}
		CSeqs Seqs;
		// Load Yeast protein sequences
		Seqs.LoadFastaFile(pro);
		// Load current positives
		Seqs.LoadPositives(pos);
		// Discover motifs in positives
		Seqs.FindingMotifs(mot);
	}
	return(0);
}
