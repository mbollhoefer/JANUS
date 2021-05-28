function [x,flag,iter,resvec]=janussolver(A,b,restart,tol,maxit,PREC,x0)
% [x,flag,iter,resvec]=janussolver(A,b,restart,tol,maxit,PREC,x0)
% [x,flag,iter,resvec]=janussolver(A,b,restart,tol,maxit,PREC)
%
% solver: Krylov subspace solver preconditioned by block ILU as computed by JANUS
%   1) Hermitian positive definite case CG
%   2) Hermitian/symmetric indefinite cas SQMR
%   3) genertal case: restarted GMRES
%
%   WARNING! make sure that A and PREC have the same structure (real/nonreal, 
%            Hermitian/symmetric/nonsymmetric, definite/non-definite)
%            For mixed cases use manually your most favoured method.
%
%   block ILU preconditioning, Q^T*SL*A*SL*P~BL*inv(BiD)*BUT^T.
%   x = janussolver(A,b,restart,tol,maxit,PREC,x0) attempts to solve
%   the system of linear equations A*x=b for x. The n-by-n coefficient matrix A
%   must be square and the right hand side column vector b must have length n.
%
%   tol specifies the tolerance of the method. 
% 
%   maxit specifies the maximum number of iterations.
% 
%   x0 optionally specifies the initial guess. If x0 is not specified then janussolver
%   uses the zero vector.
% 
%   [x,flag] = janussolver(A,b,...) also returns a convergence flag:
%     0 janussolver converged to the desired tolerance tol within maxit iterations.
%       As stopping criterion 
%       1) (CG case) the relative error in the energy norm ||x-x_k||_A/||x-x0||_A,
%       2)+3) (SQMR/GMRES case) the backward error ||Ax-b||/(||A|| ||x||+||b||)
%       is used
%     1 janussolver iterated maxit times but did not converge.
%     2 break down
%    -m unknown error code
% 
%   [x,flag,iter] = janussolver(A,b,...) also returns the iteration number
%   at which x was computed: 0 <= iter <= maxit.
% 
%   [x,flag,iter,resvec] = janussolver(A,b,...) also returns a vector of the
%   estimated residual norms at each iteration, including norm(b-A*x0).


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
% $Id: janussolver.m 4170 2018-04-26 16:35:16Z bolle $

if restart<1
   restart=30;
end

if ~isreal(A)
   PREC.isreal=0;
end

n=size(A,1);
if nargin<7
   x0=zeros(n,1);
end

if isfield(PREC,'pivots')
   pivots=PREC.pivots;
else
   pivots=[];
end   
 
% Cholesky does not need pivoting 
if PREC.isdefinite && PREC.ishermitian & PREC.invert_blocks==0
   % dummy pivot vector
   if ~isfield(PREC,'pivots')
      pivots=-(PREC.n+1);
   end
end

if PREC.isreal
   % real symmetric case
   if PREC.issymmetric
      if PREC.isdefinite
	 % hopefully the matrix is symmetric AND POSITIVE DEFINITE
	 if norm(A-A',1)==0
	    if ~isreal(b) || ~isreal(x0)
	       [xr,flag,iter,resvec]=DSPDcgbilu(A,real(b),tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,real(x0));
	       [xi,flag,iter,resvec]=DSPDcgbilu(A,imag(b),tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,imag(x0));
	       x=xr+sqrt(-1)*xi;
	    else
	       [x,flag,iter,resvec]=DSPDcgbilu(A,b,tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,x0);
	    end
	 else % nonsymmetric A
	    if ~isreal(b) || ~isreal(x0)
	       [xr,flag,iter,resvec]=Dgmresbilu(A,real(b),restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BL,PREC.P,PREC.P,PREC.SL,PREC.SL,pivots,real(x0));
	       [xi,flag,iter,resvec]=Dgmresbilu(A,imag(b),restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BL,PREC.P,PREC.P,PREC.SL,PREC.SL,pivots,imag(x0));
	       x=xr+sqrt(-1)*xi;
	    else
	       [x,flag,iter,resvec]=Dgmresbilu(A,b,restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BL,PREC.P,PREC.P,PREC.SL,PREC.SL,pivots,x0);
	    end
	 end
      else % symmetric and indefinite case
	 if norm(A-A',1)==0
	    if ~isreal(b) || ~isreal(x0)
	       [xr,flag,iter,resvec]=DSYMqmrbilu(A,real(b),tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,real(x0));
	       [xi,flag,iter,resvec]=DSYMqmrbilu(A,imag(b),tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,imag(x0));
	       x=xr+sqrt(-1)*xi;
	    else
	       [x,flag,iter,resvec]=DSYMqmrbilu(A,b,tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,x0);
	    end
	 else % non-symmetric A
	    if ~isreal(b) || ~isreal(x0)
	       [xr,flag,iter,resvec]=Dgmresbilu(A,real(b),restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BL,PREC.P,PREC.P,PREC.SL,PREC.SL,pivots,real(x0));
	       [xi,flag,iter,resvec]=Dgmresbilu(A,imag(b),restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BL,PREC.P,PREC.P,PREC.SL,PREC.SL,pivots,imag(x0));
	       x=xr+sqrt(-1)*xi;
	    else
	       [x,flag,iter,resvec]=Dgmresbilu(A,b,restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BL,PREC.P,PREC.P,PREC.SL,PREC.SL,pivots,x0);
	    end
	 end
      end % if
   % general real case
   else
      [x,flag,iter,resvec]=Dgmresbilu(A,b,restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BUT,PREC.P,PREC.Q,PREC.SL,PREC.SR,pivots,x0);
   end
else % non-real case
   % non-real complex-symmetric case
   if PREC.issymmetric
      if norm(A-A.',1)==0
	 if isreal(A)
	    [x,flag,iter,resvec]=ZSYMqmrbilu(A+(sqrt(-1)*realmin)*speye(n),b,tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,x0);
	 else
	    [x,flag,iter,resvec]=ZSYMqmrbilu(A,b,tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,x0);
	 end
      else % non-symmetric A
	 if isreal(A)
	    [x,flag,iter,resvec]=Zgmresbilu(A+(sqrt(-1)*realmin)*speye(n),b,restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BL,PREC.P,PREC.P,PREC.SL,PREC.SL,pivots,x0);
	 else
	    [x,flag,iter,resvec]=Zgmresbilu(A,b,restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BL,PREC.P,PREC.P,PREC.SL,PREC.SL,pivots,x0);
	 end % if-else	 
      end % if-else	 
   % non-real hermitian case
   elseif PREC.ishermitian
      if PREC.isdefinite
	 if isreal(A)
	    A(1,2)=A(1,2)+sqrt(-1)*realmin;
	    A(2,1)=A(2,1)-sqrt(-1)*realmin;
	    [x,flag,iter,resvec]=ZHPDcgbilu(A,b,tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,x0);
	 else
	    [x,flag,iter,resvec]=ZHPDcgbilu(A,b,tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,x0);
	 end % if-else
      else
	 if isreal(A)
	    A(1,2)=A(1,2)+sqrt(-1)*realmin;
	    A(2,1)=A(2,1)-sqrt(-1)*realmin;
	    [x,flag,iter,resvec]=ZHERqmrbilu(A,b,tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,x0);
	 else
	    [x,flag,iter,resvec]=ZHERqmrbilu(A,b,tol,maxit,PREC.BL,PREC.BiD,PREC.P,PREC.SL,pivots,x0);
	 end % if-else
      end % if-else
   % general non-real case
   else
      if isreal(A)
	 [x,flag,iter,resvec]=Zgmresbilu(A+(sqrt(-1)*realmin)*speye(n),b,restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BUT,PREC.P,PREC.Q,PREC.SL,PREC.SR,pivots,x0);
      else
	 [x,flag,iter,resvec]=Zgmresbilu(A,b,restart,tol,maxit,PREC.BL,PREC.BiD,PREC.BUT,PREC.P,PREC.Q,PREC.SL,PREC.SR,pivots,x0);
      end
   end
end
