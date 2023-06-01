#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "motif.h"

int main() {
    int length=30;

    int p0 = pow(2,0); /* first position set to 1 */
    int p1 = pow(2,2) + pow(2,1); /* second position set to 1 */
    unsigned long long int plast = pow(2,length);  /* last position set to 1 */
    unsigned long long int decimal = p0 + p1 + plast;



    printf("motif length %d\n",length);
    printf("2^%d = %lld\n", length,decimal);
    char* binaryString = ToBinary(decimal,length);
    char* Nt = ToBinary(p0,length);
    char* Nt1 = ToBinary(p1,length);
    char* Ct = ToBinary(plast,length);
    
    


    if (binaryString != NULL) {
        printf("Nter bit representation (digits=%d): %s\n", (int) strlen(Nt), Nt) ;
        printf("position (Nter+1) bit (digits=%d): %s\n", (int) strlen(Nt1), Nt1) ;
        printf("Cter representation (digits=%d): %s\n", (int) strlen(Ct), Ct) ;


        printf("Binary representation (digits=%d): %s\n", (int) strlen(binaryString), binaryString) ;
        free(binaryString);
    }


    return 0;
}