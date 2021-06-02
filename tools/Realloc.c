#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <janus.h>

#include <ilupackmacros.h>

void *REALLOC(void *ptr, size_t nbytes, char *msg )
{ /* allocates space -- or exits in case of
     failure */

  // nothing to do
  if (nbytes==0 && ptr==NULL)
     return NULL;

  // the pointer does not exist yet, do a normal allocate
  if (ptr==NULL)
     ptr=(void *)malloc(nbytes);
  // the pointer already exists and thus has to map to some existing memory 
  // region, we require a maybe different nonzero chunk of memory
  else if (nbytes!=0) // && ptr!=0
     ptr=(void *)realloc(ptr,nbytes);
  // the pointer exists and maps to some existing memory region, but we 
  // require 0 elements now. This means that we return everything and use the 
  // NULL pointer for this reason
  else //nbytes==0 && ptr!=0
     ptr=(void *)FREE(ptr);

  // we did not get the desired memory
  if (ptr==NULL && nbytes!=0) {
     fprintf(stderr,"Mem. alloc. ERROR in %s. Requested bytes: %ld bytes\n",
	     msg, (long)nbytes );
     fflush(stderr);
     exit( -1 );
  } // end if
  return ptr;
}
