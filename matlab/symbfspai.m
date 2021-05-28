function iA=symbfspai(P,tol)
% iA=symbfspai(P,tol)
% compute selected inverse based block incomplete BL BiD^{-1} BL^* factorization
% as provided by JANUS
%  
% Input:
% -----  
% P         block incomplete BL BiD^{-1} BL^* factorization WITH INVERTED DIAGONAL
%           BLOCKS (see janus for details)
% tol       accuracy of the approximate inverse, recommended to be at least as large
%           as the drop tolerance for the block ILU
%
% Output:
% -------  
% iA        selected inverse with same sparsity pattern as P(BL+BiD+BL^T)P^T, where
%           P refers to the reordering initially applied when computing
%           P^T  S_L A S_L  P ~ BL BiD^{-1} BL^* 

% $Id$

if ~P.issymmetric && ~P.ishermitian
   error('P must be symmetric/hermitian')
end

if ~P.invert_blocks
   error('turn on "invert_blocks" options to ensure that the block diagonal factor is inverted')
end

if P.isreal
  iA=DSYMbfspai(P,tol);
  n=size(iA,1);
  iA=(iA-spdiags(diag(iA),0,n,n))+iA';
elseif P.issymmetric
   iA=ZSYMbfspai(P,tol);
   n=size(iA,1);
   iA=(iA-spdiags(diag(iA),0,n,n))+iA.';
else
   iA=ZHERbfspai(P,tol);
   n=size(iA,1);
   iA=(iA-spdiags(real(diag(iA)),0,n,n))+iA';
end
