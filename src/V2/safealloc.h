/* safealloc.h */
#ifndef SAFEALLOC_H
#define SAFEALLOC_H

/* PROTOTYPES */
void*  safemalloc(size_t);         /* Allocate memory without initialization of values */
void*  safecalloc(size_t, size_t); /* Allocate and initialize memory to 0              */
void*  saferealloc(void*,size_t);  /* Reallocate memory                                */
void   safefree(void**);           /* Deallocate memory safely (set to null pointer)   */

#endif