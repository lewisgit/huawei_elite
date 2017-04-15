function [obj,x,y] = revised_simplex(c,A,b,eps,B)

% Solves: Maximize c^Tx subject to Ax = b, x >= 0
% We will assume that the LP is nondegenerate
% We are given an initial feasible basis B
% [obj,x,y] = revised_simplex(c,A,b,eps,B)
% eps is a suitable optimality tolerance, say 1e-3
% The method implements the revised simplex method in Box 7.1 on page 103
% of Chvatal
% Output parameters:- obj is the optimal objective value
% x is the primal optimal solution
% y is the dual optimal solution

% Kartik's MATLAB code for MA/OR/IE 505
% February 16, 2006
% Last modified by Kartik on April 13, 2006

[m,n] = size(A);

%%%%%%%%%%%%%%%%%%%
% Step 1:- We are given an initial basis B

N = setdiff([1:n],B);
% B = find(x0);
% N = find(ones(n,1) - abs(sign(x0)));
xB = A(:,B)\b;
% xB = x0(B);

iter = 0;
while 1 == 1,
iter = iter + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2 :- Solve A_B^Ty = c_B and compute s_N = c_N - A_N^Ty
% Declare optimality if s_N <= 0
% Else find the entering non-basic variable x_{N(k)}

y = A(:,B)'\c(B);
sN = c(N) - A(:,N)'*y;

[sNmax,k] = max(sN);

if sNmax <= eps,
   fprintf('We are done\n');
   fprintf('Number of iterations is %d\n',iter);
   x = zeros(n,1);
   x(B) = xB;
   fprintf('Optimal objective value is %f\n',c'*x);
   obj = c'*x;
   return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3 :- Solve A_Bd = a_{N(k)}
% Find theta = Min_{i=1,...,m|d_i > 0} xB(i)/d(i)
% Let theta = xB(l)/d(l)
% x_{B(l)} is the leaving basic variable
% Also check for unboundedness if d <= 0

d = A(:,B)\A(:,N(k));
zz = find(d > eps)';
if (isempty(zz))
  error('System is unbounded\n');
end
  
[theta,ii] = min(xB(zz)./d(zz));
l= zz(ii(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4:- Update B and N
% Also x(B(i)) = x(B(i)) - theta*d(i), i=1,...,m and i not equal to l
% x(B(l)) = theta

temp = B(l);
B(l) = N(k);
N(k) = temp;
xB = xB - theta*d;
xB(l) = theta;

end; % while 




