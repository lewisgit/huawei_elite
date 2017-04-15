function [x,fval,stat] = simplex(inA,inb,inp,inlb,inub)
% syntax: x = simplex(A,b,p,lb,ub)
% A two-phase simplex routine for min p'x st Ax<=b, lb<=x<=ub.
stat=1;
x=inlb;
fval=0;
A=inA;
b=inb;
p=inp;
lb=inlb;
ub=inub;
% addsize=size(A,1);
% A=[A,eye(addsize)];
% p=[p;zeros(addsize,1)];
% lb=[lb;zeros(addsize,1)];
% ub=[ub;inf*ones(addsize,1)];


if nargin < 3, error('require three input arguments'); end;
if nargin < 5, ub = inf*ones(size(p));
  if nargin < 4, lb = zeros(size(p)); end; end;
[m,n] = size(A); zer_tol = 1.0e-5; piv_tol = 1.0e-8;
 costF=p;
% set up Phase I problem

% A = sparse(A); 


v = zeros(n,1);
lbdd = find(lb>-Inf); v(lbdd) = lb(lbdd); free = setdiff(1:n,lbdd);
ubdd = free(find(ub(free)<Inf)); v(ubdd) = ub(ubdd); 
if isempty(ubdd), 
    N = -lbdd'; 
else
  N = [-lbdd' ubdd]; free = setdiff(free,ubdd);
end;
f = length(free);
% if f>0
% % %   [L,U,P,Q] = lu(A(:,free));
% % %   free = free*Q';
% % %   i = length(find(abs(diag(U))>piv_tol)); j = free(i+1:f);
% % %   if ~isempty(j)
% % %     N = [N j]; free = setdiff(free,j); f = length(free);
% % %     lb(j) = zeros(size(j)); ub(j) = zeros(size(j));
% % %   end;
% % %   [k,i] = find(P(1:f,:)); % relies on diag(U) having zeros at end
% % %   v(free) = U(1:f,1:f)\(L(1:f,1:f)\(P(1:f,:)*(b - A*v)));
% % %   k = setdiff(1:m,k); m = m - f; B = [free n+1:n+m];
% % %   
% % %   
% % %   tmp=sign(b(k)-A(k,:)*v+eps);
% % %   tmpMat=zeros(m+f,m);
% % %   for index=1:m
% % %     tmpMat(k(index),index)=tmp(index);
% % %   end
% % %   A=[A tmpMat];
%   A = [A sparse(k,1:m,sign(b(k)-A(k,:)*v+eps),m+f,m)];
  
% else
   
  B = n+1:n+m; j = [];
  tmp2=sign(b-A*v+eps*ones(size(b)));
   tmpMat2=zeros(m,m);
   for index=1:m
       tmpMat2(index,index)=tmp2(index);
   end
   A=[A tmpMat2];
%   A = [A sparse(1:m,1:m,sign(b-A*v+eps*ones(size(b))))];
% end;
lb = [lb; zeros(m,1)]; ub = [ub; Inf*ones(m,1)]; w = [zeros(n,1); ones(m,1)]; 
[x,B,N,stat] = rsmbdd(A,b,w,lb,ub,B,N);
if stat==-1
    return;
end
if (w'*x > zer_tol)
    stat=-1;
    return;
%     error('problem is infeasible'); 
end;

% convert to Phase II problem
ub(n+1:n+m) = zeros(m,1); p = [p; zeros(m,1)];
if ~isempty(j) 
  c = p(j)'-p(B)'*inv(A(:,B))*A(:,j); % should use LU decomp from rsmbdd
  if (norm(c,inf) > zer_tol)
      stat=-1;
      return;
%       error('problem is unbounded'); 
  end;
end;

[x,B,N,stat] = rsmbdd(A,b,p,lb,ub,B,N); x = x(1:n); 
if stat==-1
    return;
end
fval=costF'*x;
return;

