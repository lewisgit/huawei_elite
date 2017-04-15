% [x,val,status]=MILP_bb(c,A,b,M,eps)
% this function solves a mixed-integer linear programming problem
% using the branch and bound algorithm
% The code uses Kartik's revised simplex and dual simplex codes
% to solve the LP relaxations at each node of the branch and bound tree.
%  max c'*x
%  subject to
%        A*x = b
%        x >= 0
%        M is a vector of indices for the variables that are constrained to be integers
%        eps is the user tolerance
% the output variables are :
% x : the optimal solution
% val: value of the objective function at the optimal solution
% status =1 if successful
%        =2 if the original problem is infeasible
% Kartik's MATLAB code
% Last updated: 19th April, 2006
% The B&B driver routine based on a MATLAB code by Sherif Tawfik, Cairo University
% downloaded from the MATLAB File Exchange.

function [x,val,status]=MILP_bb(c,A,b,M,eps)

[m,n] = size(A);
% Use the phase 1 procedure to find an initial basic feasible solution
% for the LP relaxation at the root node.
Aphase1 = A;
bphase1 = b;
index = find(b < 0);
bphase1(index) = -b(index);
Aphase1(index,:) = -A(index,:);
Aphase1 = [Aphase1 eye(m)];
cphase1 = [zeros(n,1); -ones(m,1)];
x0phase1 = [zeros(n,1); bphase1];

B0 = [n+1:n+m]';

[val,x0,B] = revised_simplex_bb(cphase1,Aphase1,bphase1,eps,x0phase1,B0);

clear Aphase1 cphase1 bphase1 x0phase1 B0 index;

% If the root LP relaxation is infeasible, then just quit!
if val < 0,
  x = [];
  status = 2;
  val = [];
  return;
end  

bound = -inf; % the initial lower bound is set to -ve infinity

root = 1;

[x,B,status,b]=branch(c,A,b,x0,B,M,eps,bound,root); % a recursive function that processes the BB tree
 
val = c'*x(1:n); % objective value of the optimal solution

function [xx,B,status,bb]=branch(c,A,b,x,B,M,eps,bound,root) 
% x is an initial solution

% Solve the LP relaxation at the current node
% If it is the LP relaxation at the root node use the revised simplex method
% Else use the dual simplex method
if root == 1,
  [val0,x0,B] = revised_simplex_bb(c,A,b,eps,x,B);
  % The root LP relaxation is solved using Kartik's revised simplex code
  % Initial bfs is the one found using the phase 1 procedure.
  status0 = 1;
  root = 0; % We are now dealing with descendents of the root node
else
  [val0,x0,B,status0] = dual_simplex_bb(c,A,b,eps,x,B);
  % All LP relaxations at descendent nodes are solved using Kartik's dual
  % simplex code with warm-start information.
end   
  
% If the LP relaxation is infeasible then PRUNE THE NODE BY INFEASIBILITY
% If the LP objective value is less than our current integer bound then PRUNE THE NODE BY BOUNDS

if status0 == 2 | val0 < bound  
    xx=x; 
    status=status0; bb=bound;
    % We are pruning this node, so there is no need to branch further on
    % this node
    return;
end

% If the solution to the LP relaxation is feasible in the MILP problem, then check its objective value 
% against that of the incumbent solution
% If the new feasible solution has a greater objective value then update the lower bound
% In any case, PRUNE THE NODE BY OPTIMALITY

ind=find(abs(x0(M)-round(x0(M))) > eps ); 
if isempty(ind)
    status=1;        
    if val0 > bound % replace the incumbent solution
        x0(M)=round(x0(M));
        xx=x0;       
        bb=val0;
    else
        xx=x;  % we have pruned the node and there is no need to branch further on this node
        bb=bound;
    end
    return;
end

% The solution of the LP relaxation is not feasible in the MILP problem.
% However, the objective value of the LP relaxation is greater than the current lower bound.
% So we branch on this node to create two subproblems.
% We will solve the two subproblems recursively by calling the same branching function.
 
% The first LP problem with the added constraint that x_i <= floor(x_i) , i=ind(1)
br_var=M(ind(1));
br_value=x0(br_var);

[m,n] = size(A);
A1 = [A zeros(m,1); zeros(1,n+1)];
A1(m+1,br_var) = 1;
A1(m+1,n+1) = 1;
b1 = [b; floor(br_value)];
c1 = [c; 0];
x01 = [x0; b1(m+1)-x0(br_var)];
% The basis size grows by 1 with the new slack variable in the basis
B01 = [B; n+1];

% second LP problem with the added constraint that x_i >= ceil(x_i) , i=ind(1)
A2 = [A zeros(m,1); zeros(1,n+1)];
A2(m+1,br_var) = 1;
A2(m+1,n+1) = -1;
b2 = [b; ceil(br_value)];
c2 = [c; 0];
x02 = [x0; x0(br_var)-b2(m+1)];
% The basis size grows by 1 with the new slack variable in the basis
B02 = [B; n+1];

% solve the first child problem
[x1,B1,status1,bound1]=branch(c1,A1,b1,x01,B01,M,eps,bound,root);
status=status1;
if status1 == 1 & bound1 > bound % if the solution was successful and gives a better bound
   xx=x1;
   B=B1;
   bound=bound1;
   bb=bound1;
else
   xx=x0;
   bb=bound;
end
    
% solve the second child problem
[x2,B2,status2,bound2]=branch(c2,A2,b2,x02,B02,M,eps,bound,root);

if status2 == 1 & bound2 > bound % if the solution was successful and gives a better bound
  status=status2;
  xx=x2;
  B=B2;
  bb=bound2;
end
