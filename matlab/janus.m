function [PREC,options]=janus(A,options)
% PREC=janus(A)
% compute block incomplete L D U factorization using the
% default options
%
% This routine approximately computes Q^T SL A SR P ~  L D U,
% where SL, SR are real diagonal (scaling) matrices, P,Q are permutation 
% matrices, L is lower block unit triangular matrix and D is a block diagonal
% matrix and U upper block unit triangular
% For efficiency reasons, BL, BiD, BUT refer to cell arrays containing the 
% block structures, moreover we return BiD=D^{-1} and BUT=U^T
%
% PREC=janus(A,options)
% same routine using your private options
%
% [PREC,options]=janus(A,options)
% same routine returning the options used and possible warnings/error codes
%
%
% Input
% -----
% A        nonsingular matrix
% options  optional options. If not passed, default options are used
%          0) options.isdefinite             tell JANUS that your matrix is
%                                            (Hermitian and) positive definite
%          1) options.matching:              use maximum weight matching prior
%                                            to any reordering or factorization
%                                            default 1 (i.e., turned on),
%                                            ignored if the matrix is Hermitian
%                                            and positive definite 
%          2) options.ordering               reorder the system according to
%                                            some symmetric reordering. This
%                                            is performed preserving the structure
%                                            obtained by maximum weight matching
%                                            'mtmetis' nested dissection by nodes
%                                                      multi-threaded (default)
%                                            'metisn' nested dissection by nodes
%                                            'metise' nested dissection by edges
%                                            'amd' approximate mininum degree
%                                            'rcm' reverse Cuthill-McKee
%                                            'none' (only available if matching is
%                                                    turned off)
%          3) options.droptol                 threshold for ILU, default 1e-3
%          4) options.cosine                  use cosine preprocessing for the
%                                             initial matrix to detect dense 
%                                             blocks in the initial matrix.
%                                             default 0 (turned off)
%          5) options.blocking_strategy       'none' no further blocking is
%                                               done after cosine, matching and
%                                               reordering 
%                                             'ilu1t' default, use ILU(1,droptol)
%                                               prior to the block incomplete 
%                                               factorization to detect potential
%                                               dense blocks produced during the
%                                               incomplete factorization.
%                                             'ilupt', use ILU(p,droptol)
%                                               prior to the block incomplete 
%                                               factorization to detect potential
%                                               dense blocks produced during the
%                                               incomplete factorization, where
%                                               p refers to the level of fill
%                                             'supernodes' use supernodal
%                                               approach from direct methods to
%                                               detect dense blocks that would
%                                               occur if droptol=0
%                                             'appsuper', approximate supernodes 
%          6) options.progressive_aggregation use progressive aggregation during
%                                             the incomplete factorization to 
%                                             enlarge already existing blocks.
%                                             default 1 (turned on)
%          7) options.blocksize               optional vector predefining a 
%                                             user-defined block structure for
%                                             the diagonal blocks
%          8) options.perturbation            allow perturbation of diagonal
%                                             blocks. In the SPD case this refers
%                                             to the Ajiz&Jennigs strategy,
%                                             otherwise diagonal blocks may be
%                                             slightly perturbed by the order
%                                             of the drop tolerance
%                                             default 1 (turned on)
%          9) options.symmetric_structure     only use symmetric permutations
%                                             and attempt to incorporate 2x2 
%                                             blocks in the blocking strategy
%                                             default 0 (turned off)
%         10) options.level_of_fill           level of fill, only used when 
%                                             'ilupt' is chosen as blocking
%                                             strategy
%                                             
% 
%
% Output
% -----
% PREC       preconditioner structure
%
%  PREC.BL   block lower triangular matrix with unit diagonal stored as
%            cell array according to the number of blocks with indices
%            I (rows), J (columns), diagonal block D, sub-diagonal block E
%
%  PREC.BiD  block diagonal matrix with the same data structure as BL. The
%            diagonal blocks refer to the inverse diagonal block of D in LDU
%
%  PREC.BUT  block lower triangular matrix with unit diagonal with the same 
%            data structure as BL. BUT refers to U^T, empty in the symmetric
%            cases 
%
%  PREC.P    permutation matrices, left and right, Q empty in the symmetric
%  PREC.Q    cases
%
%  PREC.SL   real diagonal left/right scaling matrices, PREC.SR empty for
%  PREC.SR   the symmetric cases
%
%  PREC.logdet
%
% options  optional output, see input


%    Authors:
%
%	Matthias Bollhoefer, TU Braunschweig
%
%    Date:
%
%	August 28, 2017. JANUS Block ILU R1.0.  
%
%    Notice:
%
%	Copyright (c) 2017 by TU Braunschweig.  All Rights Reserved.
%
%	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
%	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
%
%    Availability:
%
%	This file is located at
%
%	http://bilu.tu-bs.de/
% $Id: janus.m 6109 2020-04-01 07:12:23Z bolle $

myoptions.isdefinite=0;
myoptions.perturbation=1;
myoptions.symmetric_structure=0;
myoptions.matching=1;
myoptions.ordering='mtmetis';
myoptions.droptol=1e-3;
myoptions.cosine='no';
myoptions.blocking_strategy='ilupt';
myoptions.level_of_fill=3;
myoptions.progressive_aggregation='yes';
myoptions.invert_blocks=1;

if nargin==2
   if isfield(options,'isdefinite')
      if isstr(options.isdefinite)
	 if strncmpi('y',options.isdefinite,1)
	    myoptions.isdefinite=1;
	 else
	    myoptions.isdefinite=0;
	 end
      elseif options.isdefinite
	 myoptions.isdefinite=1;
      else
	 myoptions.isdefinite=0;
      end
   end
   if isfield(options,'perturbation')
      if isstr(options.perturbation)
	 if strncmpi('y',options.perturbation,1)
	    myoptions.perturbation=1;
	 else
	    myoptions.perturbation=0;
	 end
      elseif options.perturbation
	 myoptions.perturbation=options.perturbation;
      else
	 myoptions.perturbation=0;
      end
   end
   if isfield(options,'symmetric_structure')
      if isstr(options.symmetric_structure)
	 if strncmpi('y',options.symmetric_structure,1)
	    myoptions.symmetric_structure=1;
	 else
	    myoptions.symmetric_structure=0;
	 end
      elseif options.symmetric_structure
	 myoptions.symmetric_structure=1;
      else
	 myoptions.symmetric_structure=0;
      end
   end
   if isfield(options,'matching')
      if isstr(options.matching)
	 if strncmpi('y',options.matching,1)
	    myoptions.matching=1;
	 else
	    myoptions.matching=0;
	 end
      elseif options.matching
	 myoptions.matching=1;
      else
	 myoptions.matching=0;
      end
   end
   if isfield(options,'ordering')
      myoptions.ordering=options.ordering;
      if (   ~strcmp(options.ordering,'mtmetis') ...
	  && ~strcmp(options.ordering,'metisn') ...
	  && ~strcmp(options.ordering,'metise') ...
          && ~strcmp(options.ordering,'amd') ...
	  && ~strcmp(options.ordering,'rcm') ...
	  && ~strcmp(options.ordering,'none'))
	 myoptions.ordering='mtmetis';
      end     
   end
   if isfield(options,'droptol')
      myoptions.droptol =options.droptol;
   end
   if isfield(options,'blocksize')
      myoptions.blocksize=options.blocksize;
   end
   if isfield(options,'cosine')
      if isstr(options.cosine)
	 if strncmpi('y',options.cosine,1)
	    myoptions.cosine='yes';
	 else
	    myoptions.cosine='no';
	 end
      elseif options.cosine
	 myoptions.cosine='yes';
      else
	 myoptions.cosine='no';
      end
   end
   if isfield(options,'blocking_strategy')
      if strncmpi('ilu1t',options.blocking_strategy,5)
	 myoptions.blocking_strategy='ilu1t';
      elseif strncmpi('ilupt',options.blocking_strategy,5)
	 myoptions.blocking_strategy='ilupt';
	 myoptions.level_of_fill=3;
      elseif strncmpi('super',options.blocking_strategy,5)
	 myoptions.blocking_strategy='supernodes';
      elseif strncmpi('appsuper',options.blocking_strategy,8)
	 myoptions.blocking_strategy='appsuper';
      elseif strncmpi('no',options.blocking_strategy,2)
	 myoptions.blocking_strategy='none';
      else
	 myoptions.blocking_strategy=options.blocking_strategy;
      end % if-elseif-else
   end
   if isfield(options,'progressive_aggregation')
      if isstr(options.progressive_aggregation)
	 if strncmpi('y',options.progressive_aggregation,1)
	    myoptions.progressive_aggregation='yes';
	 else
	    myoptions.progressive_aggregation='no';
	 end
      elseif options.progressive_aggregation
	 myoptions.progressive_aggregation='yes';
      else
	 myoptions.progressive_aggregation='no';
      end
   end
   if isfield(options,'invert_blocks')
      if isstr(options.invert_blocks)
	 if strncmpi('y',options.invert_blocks,1)
	    myoptions.invert_blocks=1;
	 else
	    myoptions.invert_blocks=0;
	 end
      elseif options.invert_blocks
	 myoptions.invert_blocks=1;
      else
	 myoptions.invert_blocks=0;
      end
   end
   if isfield(options,'level_of_fill')
      myoptions.level_of_fill=options.level_of_fill;
   end
end

if isreal(A)
   % real symmetric case
   if norm(A-A',1)==0
      if myoptions.isdefinite
	 [PREC.BL,PREC.BiD,PREC.P,PREC.SL,PREC.logdet,myoptions]=DSPDbilu(A,myoptions);
	 myoptions.isdefinite=1;
	 PREC.isdefinite=1;
      else
	 [PREC.BL,PREC.BiD,PREC.P,PREC.SL,PREC.logdet,pivots,myoptions]=DSYMbilu(A,myoptions);
	 if ~myoptions.invert_blocks
	    PREC.pivots=pivots;
	 end
	 PREC.isdefinite=myoptions.isdefinite;
      end
      PREC.issymmetric=1;
      PREC.ishermitian=1;
   else
      [PREC.BL,PREC.BiD,PREC.BUT,PREC.P,PREC.Q,PREC.SL,PREC.SR,PREC.logdet,pivots,myoptions]=DGNLbilu(A,myoptions);
      if ~myoptions.invert_blocks
	 PREC.pivots=pivots;
      end
      PREC.isdefinite=0;
      PREC.issymmetric=0;
      PREC.ishermitian=0;
   end
   PREC.isreal=1;
else 
   % non-real but complex-symmetric case
   if norm(A-A.',1)==0
      [PREC.BL,PREC.BiD,PREC.P,PREC.SL,PREC.logdet,pivots,myoptions]=ZSYMbilu(A,myoptions);
      if ~myoptions.invert_blocks
	 PREC.pivots=pivots;
      end
      PREC.isdefinite=0;
      PREC.issymmetric=1;
      PREC.ishermitian=0;
    % non-real but complex-Hermitian case
    elseif norm(A-A',1)==0
      if myoptions.isdefinite
	 [PREC.BL,PREC.BiD,PREC.P,PREC.SL,PREC.logdet,myoptions]=ZHPDbilu(A,myoptions);
	 myoptions.isdefinite=1;
	 PREC.isdefinite=1;
      else
	 [PREC.BL,PREC.BiD,PREC.P,PREC.SL,PREC.logdet,pivots,myoptions]=ZHERbilu(A,myoptions);
	 if ~myoptions.invert_blocks
	    PREC.pivots=pivots;
	 end
	 PREC.isdefinite=myoptions.isdefinite;
      end
      PREC.issymmetric=0;
      PREC.ishermitian=1;
   else
      [PREC.BL,PREC.BiD,PREC.BUT,PREC.P,PREC.Q,PREC.SL,PREC.SR,PREC.logdet,pivots,myoptions]=ZGNLbilu(A,myoptions);
      if ~myoptions.invert_blocks
	 PREC.pivots=pivots;
      end
      PREC.isdefinite=0;
      PREC.issymmetric=0;
      PREC.ishermitian=0;
   end
   PREC.isreal=0;
end % if-else
PREC.isskew=0;
PREC.n=size(A,1);

if isfield(myoptions,'isdefinite')
   if isstr(myoptions.isdefinite)
      if strncmpi('y',myoptions.isdefinite,1)
	 options.isdefinite=1;
      else
	 options.isdefinite=0;
      end
   elseif myoptions.isdefinite
      options.isdefinite=1;
   else
      options.isdefinite=0;
   end
end
if isfield(myoptions,'matching')
   if isstr(myoptions.matching)
      if strncmpi('y',myoptions.matching,1)
	 options.matching=1;
      else
	 options.matching=0;
      end
   elseif myoptions.matching
      options.matching=1;
   else
      options.matching=0;
   end
end
options.ordering=myoptions.ordering;
options.droptol =myoptions.droptol;


if isfield(myoptions,'perturbation')
   if isstr(myoptions.perturbation)
      if strncmpi('y',myoptions.perturbation,1)
	 options.perturbation=1;
      else
	 options.perturbation=0;
      end
   elseif myoptions.perturbation
      options.perturbation=myoptions.perturbation;
   else
      options.perturbation=0;
   end
end
 
if isfield(myoptions,'symmetric_structure')
   if isstr(myoptions.symmetric_structure)
      if strncmpi('y',myoptions.symmetric_structure,1)
	 options.symmetric_structure=1;
      else
	 options.symmetric_structure=0;
      end
   elseif myoptions.symmetric_structure
      options.symmetric_structure=1;
   else
      options.symmetric_structure=0;
   end
end

if isfield(myoptions,'cosine')
   if isstr(myoptions.cosine)
      if strncmpi('y',myoptions.cosine,1)
	 options.cosine=1;
      else
	 options.cosine=0;
      end
   elseif myoptions.cosine
      options.cosine=1;
   else
      options.cosine=0;
   end
end

options.blocking_strategy      =myoptions.blocking_strategy;

if isfield(myoptions,'progressive_aggregation')
   if isstr(myoptions.progressive_aggregation)
      if strncmpi('y',myoptions.progressive_aggregation,1)
	 options.progressive_aggregation=1;
      else
	 options.progressive_aggregation=0;
      end
   elseif myoptions.progressive_aggregation
      options.progressive_aggregation=1;
   else
      options.progressive_aggregation=0;
   end
end

if isfield(myoptions,'invert_blocks')
   if isstr(myoptions.invert_blocks)
      if strncmpi('y',myoptions.invert_blocks,1)
	 options.invert_blocks=1;
      else
	 options.invert_blocks=0;
      end
   elseif myoptions.invert_blocks
      options.invert_blocks=1;
   else
      options.invert_blocks=0;
   end
end
PREC.invert_blocks=options.invert_blocks;

if isfield(myoptions,'blocksize')
   options.blocksize=myoptions.blocksize;
end
if isfield(myoptions,'level_of_fill')
   options.level_of_fill=myoptions.level_of_fill;
end
