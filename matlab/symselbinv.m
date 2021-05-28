function iA=symselbinv(P)
% iA=symselbinv(P)
% compute selected inverse based block incomplete BL BiD^{-1} BL^* factorization
% as provided by JANUS
%  
% Input:
% -----  
% P         block incomplete BL BiD^{-1} BL^* factorization WITH INVERTED DIAGONAL
%           BLOCKS (see janus for details)
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
   iA=DSYMselbinv(P);
   iA=(iA+iA')./2;
elseif P.issymmetric
   iA=ZSYMselbinv(P);
   iA=(iA+iA.')./2;
else
   iA=ZHERselbinv(P);
   iA=(iA+iA')./2;
end
