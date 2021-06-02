#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <janus.h>

#include <ilupackmacros.h>

/* release memory */
void *FREE(void *ptr)
{
  /* if desired some internal statistics could be provided */
  // in real life ptr is expected to be a pointer to a pointer
  // for syntical reasons this is hidden

  if (ptr!=NULL)
     free(ptr);

  return (NULL);
} // end FREE

