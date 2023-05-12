// Class.cpp : Includes the deffinition of class Seqs
//
//#include <Windows.h>
//#include <process.h>
#include <pthread.h>
#include "Utils.hpp"

// ---------------------------------------------------------------------------
// Structur : Protein Sequences
class CSeqs
{
	private:

	public:
		// Creator Destructor
		CSeqs(); 																// Constructor
	   ~CSeqs(); 																// Destructor

		// Variables
		CSeq*	Seq; 															// The sequence
		int		SeqPOS;															// Sequence number
		int		SeqPRO;															// Sequence number

		// Routines
		void	LoadFastaFile(char *Inputfile);						// Load sequences from FASTA format
		void	LoadPositives(char *Inputfile);						// Load sequences from FASTA format
		void	FindingMotifs(char *OutputFile);					// Detect all motifs in positives
//		void	TestMotifs(char *OutputFile);
};
// ---------------------------------------------------------------------------
// Constructor
CSeqs::CSeqs()
{
	SeqPOS = 0;
	SeqPRO = 0;
	Seq    = NULL;
}
// ---------------------------------------------------------------------------
// Destructor
CSeqs::~CSeqs()
{
	for (int i=0; i<SeqPRO; i++)
	delete[] Seq[i].AA;
	delete[] Seq;
}
// ---------------------------------------------------------------------------
// Load sequences from FASTA format
void CSeqs::LoadFastaFile(char *Inputfile)
{
	// Verbose
	printf("Loading proteome\n");

	int   i;
	char  Buff[MAXLEN];
	FILE *file;

	// On ouvre le fichier en lecture
	if ((file = fopen(Inputfile, "r")) == NULL)
	{
		fprintf(stderr,"# Input file not found ! #");
		return;
	}

	// Chercher le nombre de séquences
	SeqPRO = 0;
	while (!feof(file))
	{
		// read the current line from the input file
		if (fgets(Buff, MAXLEN, file) == NULL) continue;

		// Increment the number of sequences
		if (Buff[0] == '>') SeqPRO++;
	}

	// Retourner au début du fichier
	fseek(file, 0, SEEK_SET);

	// Initialiser le nombre de séquences
	Seq = new CSeq[SeqPRO];

	// Chercher la longueur de chaque séquence
	SeqPRO = 0;

	while (!feof(file))
	{
		// read the current line from the input file
		if (fgets(Buff, MAXLEN, file) == NULL) continue;

		// Clean the buffer
		strtok(Buff,"\n");
		strtok(Buff,"\r");
		strtok(Buff," ");

		if (Buff[0] == '>')
		{
			// Initialiser le vecteur de séquences
			Seq[SeqPRO].Num    = SeqPRO;
			Seq[SeqPRO].TYP    = SEQUENCE;
			Seq[SeqPRO].Length = 0;

			// Get the name of the sequence
			strcpy(Seq[SeqPRO].Name, strtok(Buff+1, "\t"));

			// Increment sequences number
			SeqPRO++;

			// Go to the next line
			continue;
		}

		// Increment the length of the sequence
		Seq[SeqPRO-1].Length += strlen(Buff);
	}

	// Intialisation des vecteur séquences et fréquences
	for (i=0; i<SeqPRO; i++)
	{
		Seq[i].AA = Alloc1D(Seq[i].Length+1, ' ');

		// Add terminate character
		Seq[i].AA[Seq[i].Length] = '\0';
	}

	// Retourner au début du fichier
	fseek(file, 0, SEEK_SET);

	// Lire le contenu de chaque séquence
	SeqPRO = 0;

	while (!feof(file))
	{
		// read the current line from the input file
		if (fgets(Buff, MAXLEN, file) == NULL) continue;

		// Remove the end line caracter
		strtok(Buff,"\n");
		strtok(Buff,"\r");
		strtok(Buff," ");

		if (Buff[0] == '>')
		{
			// Initialize the protein length
			Seq[SeqPRO++].Length = 0;
			continue;
		}

		for (i=0; i<strlen(Buff); i++)
		{
			if ((Buff[i]) == '.' || (Buff[i]) == '~') (Buff[i]) = '-';
			if ((Buff[i]) == '*' || (Buff[i]) == 'U') (Buff[i]) = 'X';
			if (isalpha(Buff[i]) || (Buff[i]) == '-') Seq[SeqPRO-1].AA[Seq[SeqPRO-1].Length++] = toupper(Buff[i]);
		}
	}
	// Fermer le fichier
	fclose(file);

	// Information sur le fichiers
	printf("Input file          : %s\n", Inputfile);
	printf("Number of sequences : %d\n", SeqPRO);
}
// ---------------------------------------------------------------------------
// Load positives sequences
void CSeqs::LoadPositives(char *Inputfile)
{
	// Verbose
	printf("Loading positives\n");

	CSeqs Temp;

	// Load Yeast protein sequences
	Temp.LoadFastaFile(Inputfile);

	// Check if positives loaded
	if (Temp.SeqPRO == 0) return;

	for (int i=0; i<Temp.SeqPRO; i++)
	for (int j=0; j<SeqPRO; j++)
	{
		if (strcmp(Temp.Seq[i].Name, Seq[j].Name) == 0)
		{
			if (Temp.Seq[i].Length != Seq[j].Length)
			{
				fprintf(stderr,"# Difference in length of sequence in positives and proteome ! #");
				exit(EXIT_FAILURE);
			}

			// Increment positives
			SeqPOS += 1;

			// Set sequence to positive
			Seq[j].TYP = POSITIVE;

			// Copy positive sequence
			strcpy(Seq[j].AA, Temp.Seq[i].AA);

			// Go next
			break;
		}
	}

	// Just a test
	if (Temp.SeqPRO != SeqPOS)
	{
		fprintf(stderr,"# Positives not in proteome ! #");
		exit(EXIT_FAILURE);
	}
}
// ---------------------------------------------------------------------------
// Discover motifs in positives
void CSeqs::FindingMotifs(char *OutputFile)
{
	// Verbose
	printf("Discover motifs\n");

	char * inutile;
	char   BUF[MAXLEN];
	FILE  *file1;
	FILE  *file2;

	// Open input and output files
	file1 = fopen(OutputFile, "w"); fclose(file1);

	for (int I=MINMOT; I<=MAXMOT; I++)
	{
		printf("%d\n",I);
		for (int J=0; J<20; J+=THREAD)
		{
			int threadId;
			pthread_t * hThreads= new pthread_t [THREAD];

			CThread **Thread = new CThread * [THREAD];
			for (int K=0; K<THREAD; K++)
			{
				Thread[K]      = new CThread[1];
				Thread[K]->LEN = I;
				Thread[K]->RES = J+K;
				Thread[K]->THR = K;
				Thread[K]->POS = SeqPOS;
				Thread[K]->PRO = SeqPRO;
				Thread[K]->SEQ = Seq;
			}

			for (int K=0; K<THREAD; K++)
			{
				// Create Thread
				threadId = pthread_create( &hThreads[K], NULL, MotifsThread, (void*) Thread[K]);

				// Check if thread successfuly created
				if(threadId != 0)
				{
					fprintf(stderr,"Error! Creation Failed for Thread %d\n",threadId);
					exit(1);
				}
			}

			// Resume suspended threads
			for (int K=0; K<THREAD; K++)
			pthread_join( hThreads[K], NULL);

			// Empty current threads
			for (int K=0; K<THREAD; K++)
			delete[] Thread[K];
			delete[] Thread;
			delete[] hThreads;

			// Open input and output files
			if ((file1 = fopen(OutputFile, "a")) == NULL)
			{
				fprintf(stderr,"# Ouput filepath not found ! #");
				exit(EXIT_FAILURE);
			}

			for (int K=0; K<THREAD; K++)
			{
				// Open current output file
				sprintf(BUF, "TMP_%d.dat", K);
				file2 = fopen(BUF, "r");

				while (!feof(file2))
				{
					// read the current line from the input file
					if (fgets(BUF, MAXLEN, file2) == NULL) continue;

					// Save current motif
					fprintf(file1, "%s", BUF);
				}

				// close current stream
				fclose(file2);
			}
			// Close output stream
			fclose(file1);
		}
	}
	// Delete temporary files
	printf("End Discover motifs\n");
	system("rm TMP*.dat");
}
/// ---------------------------------------------------------------------------
//// Test motifs
//void CSeqs::TestMotifs(char *OutputFile)
//{
//	// Verbose
//	printf("Test motifs\n");

//	int    NBR=0;
//	char   MOT[MAXLEN];
//	char   BUF[MAXLEN];
//	int    POS1;
//	int    POS2;
//	int    POS3;
//	int    PRO1;
//	int    PRO2;
//	int    PRO3;
//	double PVAL;
//	FILE *file;

//	unsigned int LEN;
//	unsigned int OFF;

//	// Open input file
//	file = fopen(OutputFile, "r");

//	while (!feof(file))
//	{
//		// Save motif statistics
//		if (fscanf(file, "%s\t%d\t%d\t%d\t%d\t%lf", &MOT, &POS1, &PRO1, &POS2, &PRO2, &PVAL) == EOF) continue;

//		// Verbose
//		printf("%d\r", ++NBR);

//		// Init occurrences
//		POS3 = 0;
//		PRO3 = 0;

//		// Prepare regular expression
//		TRegexp RegExp(MOT);

//		for (int i=0; i<SeqPRO; i++)
//		{
//			if ((OFF = RegExp.find(Seq[i].AA, &(LEN=0), (OFF=0))) <= Seq[i].Length - LEN)
//			{
//				// Occurrence in proteome
//				POS3 += Seq[i].TYP == POSITIVE;
//				PRO3 += 1;
//			}
//		}

//		// Free regexp
//		RegExp.~TRegexp();

//		// Check error
//		if (POS1 != POS3 || PRO1 != PRO3) printf("%s\n", MOT);
//	}
//	// Close stream
//	fclose(file);
//}

