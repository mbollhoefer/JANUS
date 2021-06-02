/* $Id: janusdefault.c 6291 2020-05-27 20:25:07Z bolle $ */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <blas.h>
#include <lapack.h>

#include <janus.h>
#include <ilupackmacros.h>

void JanusDefaultOptions(JanusOptions *options) {

    options->matching=1;                /* turn matching on/off */
    options->ordering=PERM_MTMETIS;     /* which reordering to take */
    options->droptol=1e-3;              /* drop tolerance  */
    options->cosine=0;                  /* cosine-based blocking turned on/off */
    options->blocking_strategy=BLOCK_ILUPT; /* ILU(p,droptol) blocking strategy turned on */
    options->progressive_aggregation=1; /* progressively aggreate blocks during the factorization on/off */
    options->perturbation=1;            /* allow diagonal block perturbations */
    options->blocksize=NULL;            /* optionally, pass a fixed block partitioning*/
    options->nblocks=0;                 /* number of blocks */
    options->symmetric_structure=0;     /* use symmetric permutations when the nonsymmetric matrix has almost symmetric structure */
    options->invert_blocks=1;           /* invert diagonal blocks */
    options->level_of_fill=5;           /* level of fill, if ILUPT is chosen */
}
