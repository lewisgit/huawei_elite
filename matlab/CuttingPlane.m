
%this function solves a mixed-integer linear programming problem
% using the gomery cutting plane algorithm
% The code uses MATLAB's linear programming solver "linprog"
% to solve the LP relaxations at each node of the branch and bound tree.
%   min f*x
%  subject to
%        A*x <=b
%        Aeq * x = beq   
%        lb <= x <= ub
%        M is a vector of indeces for the variables that are constrained to be integers
%        B is a vector of indeces for the basic variables
%        itermax is the variable for the maximum number of gomery cutting
%                   planes added
%
% 
% By Shon Albert, Department of Mathematics, North Carolina State University
% 
% Last updated: 12th November, 2006.

function [c,A,b,obj,x,B,iter] = CuttingPlane(c,A,b,eps,x0,B,M,itermax)


%
%first find the objective value of x0 
%
obj=x0*c;
%
%solve for the size of A
%
[m,n]=size(A);
p=n;
%
% Use the phase 1 procedure to find an initial basic feasible solution
% for the LP relaxation at the root node.
%
Aphase1 = A;
bphase1 = b;
index = find(b < 0);
bphase1(index) = -b(index);
Aphase1(index,:) = -A(index,:);
Aphase1 = [Aphase1 eye(m)];
cphase1 = [zeros(n,1); -ones(m,1)];
x0phase1 = [zeros(n,1); bphase1];
B0 = [n+1:n+m]';

[val,x,B] = revised_simplex_bb(cphase1,Aphase1,bphase1,eps,x0phase1,B0);
%
%reduce x to its orignal size and sort B
%
x=x(1:n);
B=sort(B);
%
%solve the LP relaxation with Revised Simplex Method
%
[val,x,B] = revised_simplex_bb(c,A,b,eps,x,B);
%
%Check if the solution is integer
%if so t=1 else t=0
%
t=isequal(x(M),floor(x(M)+eps));
B=B';
%
%set iter=0 and do the following steps until t=1 or iter=itermax
%
iter=0;
while(t==0 && iter<itermax)
%
%Solve for the Nonbasic variables
%
N = setdiff([1:n],B);
%
%Solve for the optimal dictionary
%
D = ((A(:,B))\A(:,N));
%
%set br_var = most fractional
%
[C,br_var]=max(abs(x(M)-round(x(M))))
%
%set br_value
%
br_value=floor(x(br_var)+eps)
%
%get z= the gomery cutting plane and add the cut to the system
%
z=D(br_var,:);

    A = [A zeros(m,1); zeros(1,n+1)];
    m=m+1;
    n=n+1;

    A(m,br_var) = 1;    
    A(m,n) = 1;
    A(m,N)=floor(D(br_var,:));
    b = [b; br_value];
    c = [c; 0];
    x = [x; b(m)-sum(A(m,[1:n-1])*x)]
    B = [B'; n]';
   
%
%now solve using the dual simplex method and set iter=iter+1
%

    [obj,x,B,status] = dual_simplex_bb(c,A,b,eps,x,B);      
iter=iter+1
%
%check if we have an integer solution, if so t=1
%
if(x(M)-floor(x(M)+eps) < eps)
    t=1;
end
%
%end the while loop
%
end
%
%return x to its original size
%
x=x(1:p);