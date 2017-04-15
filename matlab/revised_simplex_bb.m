function [obj,x,B] = revised_simplex_bb(c,A,b,eps,x0,B)

% Solves: Maximize c^Tx subject to Ax = b, x >= 0
% We are given an initial feasible basis B and the corresponding bfs x0
% [obj,x,y] = revised_simplex(c,A,b,eps,x0,B)
% eps is a suitable optimality tolerance, say 1e-6
% The method implements the revised simplex method in Box 7.1 on page 103
% of Chvatal
% Output parameters:- obj is the optimal objective value
% x is the primal optimal solution
% B is the optimal basis for the LP

% Kartik's MATLAB code for MA/OR/IE 505
% February 16, 2006
% Last modified by Kartik on April 19, 2006
% Designed for use in the branch and bound scheme

[m,n] = size(A);

%%%%%%%%%%%%%%%%%%%
% Step 1:- We are given an initial basis B

N = setdiff([1:n],B);
xB = x0(B);

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
   x = zeros(n,1);
   x(B) = xB;
   obj = c'*x;
   status = 1;
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
if (isempty(zz)),
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




