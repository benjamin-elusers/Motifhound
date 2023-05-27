#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "safealloc.h"

/* PROTOTYPES */ 
void*  safemalloc(size_t);          /* Allocate memory without initialization of values */
void*  safecalloc(size_t, size_t);  /* Allocate and initialize memory to 0              */
void*  saferealloc(void*,size_t);   /* Reallocate memory                                */
void   safefree(void**);            /* Deallocate memory safely (set to null pointer)   */

/* FUNCTIONS */
void* safemalloc(size_t size) {
/* malloc with exit on failed memory allocation */
    void* ptr = malloc(size);
    if (ptr == NULL) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(ENOMEM);
    }
    return ptr;
}

void* safecalloc(size_t num, size_t size) {
/* calloc with exit on failed memory allocation */
    void* ptr = calloc(num, size);
    if (ptr == NULL) {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(ENOMEM);
    }
    return ptr;
}

void* saferealloc(void *ptr, size_t s) {
/* realloc with exit on failed memory allocation */
	void *newptr = realloc(ptr,s);
	if(newptr == NULL) { 
		fprintf(stderr, "Memory reallocation failed!\n");
		exit(ENOMEM);    }
	return newptr;
}

void safefree(void** ptr){
     free(*ptr);
    *ptr = NULL;
}

