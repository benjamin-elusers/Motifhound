//---------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
//#include <regexp.h>
#define MAXLEN   1024							// length of taxon names
#define MIN2(a,b)   (a < b ? a : b)
#define NBFAST		80

#define SEQUENCE    0
#define POSITIVE    1
#define NEGATIVE    2


int FRQMOT =     3;							// Minimum number of occurence to keep motif in suffix tree
int MINMOT =     3;							// Minimum length of a motif
int MAXMOT =     8;							// Maximum length of a motif
int NWCMOT =     3;							// Minimum number of non-wildcards positions in a motif
int THREAD =     4;							// Number of simultaneous threads

//---------------------------------------------------------------------------
// Structure for Suffix Tree
typedef struct CTree
{
	int     IDX;
	char    NBR;
	char    RES;
	CTree **SON;
} CTree;
//---------------------------------------------------------------------------
// Structure for sequence
typedef struct CSeq
{
	int   Num;								// Sequence index
	int   Length;							// Length of the sequence
	int   TYP;								// Sequence type
	char  Name[MAXLEN];						// Sequence name
	char *AA;								// Amino acids
} CSeq;
//---------------------------------------------------------------------------
// Structure for thread
typedef struct CThread
{
	int   THR;
	int   LEN;
	int   RES;
	int   POS;
	int   PRO;
	CSeq *SEQ;
} CThread;
//---------------------------------------------------------------------------
// Swap two variables
template<class T>
inline void SWAP(T &a, T &b)	{T t=a; a=b; b=t;}
//---------------------------------------------------------------------------
// Resize memory allocation
template<class T>
inline T* Realloc(T *ptr, int m)
{
	if ((ptr = (T *)realloc(ptr, m * sizeof(T))) == NULL)
	{
		printf("FATAL EEROR !!");
		exit(EXIT_FAILURE);
	}
	// Return new memory allocation
	return ptr;
}
//---------------------------------------------------------------------------
// Resize memory allocation
template<class T>
inline T* Alloc1D(int m, T v)
{
	// Allocate memory
	T *ptr = new T[m];

	// Init values
	for (int i=0; i<m; i++) ptr[i] = v;

	// Return new memory allocation
	return ptr;
}
//---------------------------------------------------------------------------
// Resize memory allocation
template<class T>
inline T** Alloc2D(int m, int n, T v)
{
	// Allocate memory
	T **ptr = new T*[m];

	// Init values
	for (int i=0; i<m; i++)
	{
		// Allocate memory
		ptr[i] = new T[n];

		for (int j=0; j<n; j++) ptr[i][j] = v;
	}

	// Return new memory allocation
	return ptr;
}
//---------------------------------------------------------------------------
// Get residue index from its char
inline char Char(char AA)
{
	switch  (AA)
	{
		case  0 : return 'A';
		case  1 : return 'R';
		case  2 : return 'N';
		case  3 : return 'D';
		case  4 : return 'C';
		case  5 : return 'Q';
		case  6 : return 'E';
		case  7 : return 'G';
		case  8 : return 'H';
		case  9 : return 'I';
		case 10 : return 'L';
		case 11 : return 'K';
		case 12 : return 'M';
		case 13 : return 'F';
		case 14 : return 'P';
		case 15 : return 'S';
		case 16 : return 'T';
		case 17 : return 'W';
		case 18 : return 'Y';
		case 19 : return 'V';
		default  : exit(EXIT_FAILURE);
	}
}
// ---------------------------------------------------------------------------
// Choose and update the substitution matrix
void SaveSeq(CSeq Seq, char *OutputFile, char *opt)
{
	int   i, j;
	FILE *file;

	file = fopen(OutputFile, opt);
	fprintf(file, ">%s\n", Seq.Name);

	for (j=0; j<Seq.Length; j++)
	{
		fprintf(file, "%c", Seq.AA[j]);

		div_t x = div(j+1, NBFAST);
		if (x.rem == 0 || j == Seq.Length-1) fprintf(file, "\n");
	}
	fprintf(file, "\n");
	fclose(file);
}
// ---------------------------------------------------------------------------
// Compute binomial coefficient
double BinomialCoefficient(double n, double k)
{
	double b = 1.0f;
	double x = 1.0f;

	while (x <= k)
	{
		b *= (n - (k - x)) / x;
		x += 1.0f;
	}

	return b;
}
// ---------------------------------------------------------------------------
// Compute cumulative hypergeometric distribution
double HyperGeometric(double k, double n, double m, double N)
{
	double c;
	double h;
	double x;

	c = h = BinomialCoefficient(m, k) * BinomialCoefficient(N - m, n - k) / BinomialCoefficient(N, n);

	for (x=k+1; x<=MIN2(m,n); x+=1.0f)
	{
		h *= (x - (m + 1.0f)) * (x - (n + 1.0f)) / (x * (x + (N - m - n)));
		c += h;
	}

	return c;
}
// ---------------------------------------------------------------------------
// Add a motif in suffix tree
void SetTree(CTree *Tree, char *MOT, int LEN, int POS, int NWC, int *NBR, char ***STR, char *CHR)
{
	CTree *ptr = Tree;

	// Look for a match
	while (ptr != NULL)
	{
		bool fnd = false;

		for (int i=0; i<ptr->NBR; i++)
		{
			if (ptr->SON[i]->RES == MOT[POS] || ptr->SON[i]->RES == '.')
			{
				// Check if the motif had been already added
				if (POS == LEN-1) return;

				// Get current position
				CHR[POS] = ptr->SON[i]->RES;

				// Go Next
				SetTree(ptr->SON[i], MOT, LEN, POS+1, NWC+(ptr->RES == '.'), NBR, STR, CHR);

				// Match
				if (ptr->SON[i]->RES != '.')
				{
					ptr = ptr->SON[i];
					fnd = true;
					POS++;
					break;
				}
			}
		}

		if (fnd == false)
		{
			if ((ptr->NBR == 0) && (NWCMOT < LEN-NWC) && (POS > 0) && (POS < LEN-1))
			{
				// Get current position
				CHR[POS] = '.';

				// Number of sons
				ptr->NBR = 1;

				// Create new node
				ptr->SON         = new CTree*[1];
				ptr->SON[0]      = new CTree[1];
				ptr->SON[0]->IDX = -1;
				ptr->SON[0]->NBR =  0;
				ptr->SON[0]->RES = '.';
				ptr->SON[0]->SON = NULL;

				// Go Next
				SetTree(ptr->SON[0], MOT, LEN, POS+1, NWC+1, NBR, STR, CHR);
			}

			// Get current position
			CHR[POS] = MOT[POS];

			// Number of sons
			ptr->NBR += 1;

			// Create new node
			ptr->SON                  = Realloc(ptr->SON, ptr->NBR);
			ptr->SON[ptr->NBR-1]      = new CTree[1];
			ptr->SON[ptr->NBR-1]->IDX = POS < LEN-1 ? -1 : (*NBR)++;
			ptr->SON[ptr->NBR-1]->NBR = 0;
			ptr->SON[ptr->NBR-1]->RES = CHR[POS];
			ptr->SON[ptr->NBR-1]->SON = NULL;

			if (++POS == LEN)
			{
				// Update pointers array
				(*STR) = Realloc(*STR, *NBR);
				(*STR)[(*NBR)-1] = Alloc1D(LEN+1, ' ');
				strncpy((*STR)[(*NBR)-1], CHR, LEN);
				(*STR)[(*NBR)-1][LEN] = '\0';

				// Check motif's end
				return;
			}

			// Go next
			ptr = ptr->SON[ptr->NBR-1];
		}
	}
}
// ---------------------------------------------------------------------------
// Get a motif from Suffix Tree
void GetTree(CTree *Tree, char *MOT, int LEN, int POS, int TYP, bool *FND, int **OCC)
{
	CTree *ptr = Tree;

	// Look for a match
	for (int i=0; i<ptr->NBR; i++)
	{
		if (ptr->SON[i]->RES == MOT[POS] || ptr->SON[i]->RES == '.')
		{
			// Go next
			GetTree(ptr->SON[i], MOT, LEN, POS+1, TYP, FND, OCC);

			// Check if the motif had been already added
			if (POS == LEN-1)
			{
				// Increment frequence of current motif
				if (FND[ptr->SON[i]->IDX] == false)
				{
					// Tag current motif
					FND[ptr->SON[i]->IDX] = true;

					// Get motifs occurrences
					OCC[ptr->SON[i]->IDX][0] += TYP == POSITIVE;
					OCC[ptr->SON[i]->IDX][1] += 1;
				}

				// Exit
				break;
			}

			// Skip if non wildcard match
			if (ptr->SON[i]->RES != '.') break;
		}
	}
}
// ---------------------------------------------------------------------------
// Empty memory allocated for Suffix Tree
void EmptyTree(CTree *Tree)
{
	// Empty memory
	for (int i=0; i<Tree->NBR; i++)
	EmptyTree(Tree->SON[i]);

	// Empty memory
	delete[] Tree->SON;
	delete[] Tree;
}
//// ---------------------------------------------------------------------------
//// Discover motifs in Set of proteins of interest
void * MotifsThread(void *arg)
{
	// get local pointer from arg
	CThread *ptr = (CThread*)arg;

	int    NBR;
	char   BUF[MAXLEN];
	char   CHR[MAXLEN];
	double VAL;
	FILE  *file;

	bool  *FND;		// Presence/ABcence of motif
	int  **OCC;		// Motif occurrence
	char **MOT;		// Enumerated motif

	// Prepare Suffix Tree
	CTree *Tree;

	// Amino Acid in first positoin of the current motif 
	char RES = Char(ptr->RES);

	// Number of motifs
	NBR = 0;
	MOT = Alloc2D(1, ptr->LEN, ' ');

	// Prepare Suffix Tree
	Tree      = new CTree[1];
	Tree->IDX = -1;
	Tree->NBR =  0;
	Tree->RES = '#';
	Tree->SON = NULL;

	// Build suffix tree
	for (int i=0; i<ptr->PRO; i++)
	{
		// Skip if not ligand
		if (ptr->SEQ[i].TYP != POSITIVE) continue;

		for (int j=0; j<ptr->SEQ[i].Length-ptr->LEN+1; j++)
		{
			// Check motif's first residue
			if (ptr->SEQ[i].AA[j] != RES) continue;

			// Add motif in the suffix tree
			SetTree(Tree, &ptr->SEQ[i].AA[j], ptr->LEN, 0, 0, &NBR, &MOT, CHR);
		}
	}

	// Tag motifs presence
	FND = Alloc1D(NBR, false);
	OCC = Alloc2D(NBR, 2, 0);

	// Build suffix tree
	for (int i=0; i<ptr->PRO; i++)
	{
		// Init motifs presence
		memset(FND, false, NBR * sizeof(bool));

		for (int j=0; j<ptr->SEQ[i].Length-ptr->LEN+1; j++)
		{
			// Check motif's first residue
			if (ptr->SEQ[i].AA[j] != RES) continue;

			// Add motif in the suffix tree
			GetTree(Tree, &ptr->SEQ[i].AA[j], ptr->LEN, 0, ptr->SEQ[i].TYP, FND, OCC);
		}
	}

	// Empty memory
	delete[] FND;

	// Open input and output files
	sprintf(BUF, "TMP_%d.dat", ptr->THR); file = fopen(BUF, "a");

	// Store motif's occurrences
	for (int i=0; i<NBR; i++)
	{
		// Skip if motif occurs less than FRQMOT
		if (OCC[i][0] < FRQMOT) continue;

		// Compute motif's significance
		VAL = HyperGeometric(OCC[i][0], ptr->POS, OCC[i][1], ptr->PRO);

		// Skip if not significant
		fprintf(file, "%s\t%d\t%d\t%d\t%d\t%.10e\n", MOT[i], OCC[i][0], OCC[i][1], ptr->POS, ptr->PRO, VAL);
	}

	// Close output stream
	fclose(file);

	for (int i=0; i<NBR; i++)
	{
		delete[] OCC[i];
		delete[] MOT[i];
	}
	delete[] OCC;
	delete[] MOT;

	// Empty Suffix Tree
	EmptyTree(Tree);

	// Return value
//	return 0;
}
