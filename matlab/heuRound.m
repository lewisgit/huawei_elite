function [xx,ff,stat]=heuRound(f,intcon,A,b,lb,ub)
% curit=curit+1;

[xx,ff,stat]=simplex(A,b,f,lb,ub);
frac=find(abs(xx(intcon)-0.5)<(0.5-1e-4),1);
if(isempty(frac) && stat~=-1)
    return;
end

stat=-1;
while(stat==-1)
    roundx=xx(intcon);
    roundx=roundx/max(roundx);
    roundx=roundx+0.1;
    roundxx=round(roundx);
    while(sum(roundx)>size(intcon))
        roundx=roundx-0.1;
        roundxx=round(roundx);
    end
    nlb=lb;nub=ub;
    nlb(intcon)=roundxx;
    nub(intcon)=roundxx;
    [xxx,fff,stat]=simplex(A,b,f,nlb,nub);
end
xx=xxx;
ff=fff;


