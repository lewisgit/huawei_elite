function [x,fval]=branchNbound(f,intcon,A,b,lb,ub,probsize,sol0,val0)

[x,fval] = simplex(A,b,f,lb,ub);
isInt=find(abs(x(2:probsize(1)-0.5))>1e-4);
if(isempty(isInt))
    return;
end
[~,index]=min(abs(x(1:probsize(1))-0.5));
lb1=lb;
lb1(index)=1;
ub2=ub;
ub2(index)=0;
[x1,fval1]=branchNbound(f,intcon,A,b,lb1,ub,probsize,sol0,val0);
[x2,fval2]=branchNbound(f,intcon,A,b,lb,ub2,probsize,sol0,val0);
