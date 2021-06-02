/* $Id: etree.c 4077 2018-03-16 15:37:45Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <janus.h>
#include <ilupackmacros.h>

#ifdef _USE_MKL_
#include <mymkl.h>
#else
#endif

/* compute the etree of A(p,p), adapted from Tim Davis' book

   A        sparse matrix

   p	    permutation, FORTRAN-style indexing
   invq     inverse permutation with respect to p, 
            FORTRAN-style indexing

   parent   array of parents in the elimination tree of A(p,p)

   iw       work space of length A.nc used to store ancestors
            which leads to path compression and thus more
	    efficiency

   Remark. To compute the elimination tree, the graph of the
           matrix need not be connected
           The columns of A(p,p) are required to have indices in
	   increasing order
*/
void ETREE(SPARSEMATRIX *A,
	   integer *p, integer *invq, integer *parent,
	   integer *iw)
{
    integer i,k,j,l, inext, *ancestor=iw, n=A->nc,*ja;

    for (k=0; k<n; k++) {
        // node k has no parent yet
	parent[k]=-1;		    
	// nor does k have an ancestor
	ancestor[k]=-1;
	// number of nonzero entries in column p[k]
	l=A->ncol[p[k]];
	ja=A->rowind[p[k]];
	for (j=0; j<l; j++) {
	    i=invq[ja[j]];
	    // traverse from i to k
	    while (i!=-1 && i<k) {
	          // next node = ancestor of i
	          inext=ancestor[i];		    
		  // path compression
		  ancestor[i]=k;		    
		  // no ancestor, then parent is k
		  if (inext==-1) 
		     parent[i]=k;   
		  i=inext;
	    } // end while
	} // end for j
    } // end for k
} // end ETREE






/* compute the etree of A, adapted from Tim Davis' book

   A        sparse matrix

   parent   array of parents in the elimination tree of A

   iw       work space of length A.nc used to store ancestors
            which leads to path compression and thus more
	    efficiency

   Remark. To compute the elimination tree, the graph of the
           matrix need not be connected
           The columns of A are required to have indices in
	   increasing order
*/
void ETREE0(SPARSEMATRIX *A,
	   integer *parent,
	   integer *iw)
{
    integer i,k,j,l, inext, *ancestor=iw, n=A->nc,*ja;

    for (k=0; k<n; k++) {
        // node k has no parent yet
	parent[k]=-1;		    
	// nor does k have an ancestor
	ancestor[k]=-1;
	// number of nonzero entries in column k
	l=A->ncol[k];
	ja=A->rowind[k];
	for (j=0; j<l; j++) {
	    i=ja[j];
	    // traverse from i to k
	    while (i!=-1 && i<k) {
	          // next node = ancestor of i
	          inext=ancestor[i];		    
		  // path compression
		  ancestor[i]=k;		    
		  // no ancestor, then parent is k
		  if (inext==-1) 
		     parent[i]=k;   
		  i=inext;
	    } // end while
	} // end for j
    } // end for k
} // end ETREE0

