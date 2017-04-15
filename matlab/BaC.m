%this function solves a mixed-integer linear programming problem
% using the gomory cutting plane scheme and the branch and bound algorithm
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

function [x,iter1,iter2]=BaC(c,A,b,eps,x0,B,M,itermax)
%
%if itermax = 0 then we set iter1 = 0
%else we call the CittingPlane function to make
%itermax number of gomory cuts
%
if itermax==0
    
    iter1=0;
else

[c,A,b,obj,x,B,iter1] = CuttingPlane(c,A,b,eps,x0,B,M,itermax);

end
%
%if the system is not completely solved then we solve the 
%remaining system by calling the branch and bound function
%
if(iter1==itermax)

[x,val,status,iter2]=MILP_bb(c,A,b,M,eps);


else
    
    iter2=0;
    
end

