%edges, and conReq are alreade loaded;
% probSize=[50 99 9];
% servCost=240;
mat=genMat(edges,probSize(1));
% [f,intcon,A,b,Aeq,beq,lb,ub]=genConst(probSize,servCost,edges,conReq,mat);
[f,intcon,A,b,Aeq,beq,lb,ub]=genConst(probSize,servCost,edges,conReq,mat);
%  [x,fval,exitflag,output]=linprog(f,A,b,[],[],lb,ub);
options = optimoptions('intlinprog','RootLPAlgorithm','dual-simplex','CutMaxIterations',30,'BranchRule','maxfun','HeuristicsMaxNodes',100,'MaxNodes',10000000,'Heuristics','round','IntegerPreprocess','basic','CutGeneration','basic','NodeSelection','minobj');
% [x,fval,exitflag,output]=intlinprog(f',intcon,A,b,Aeq,beq,lb,ub,options);
% [x2,fval2,stat] = simplex(A,b,f',lb,ub);
%  [xx,ff,stat2]=heuRound(f',intcon,A,b,lb,ub);
% [x,fval,status]=MILP_bb(f',A,b,intcon,1e-6);
xtype=1:probSize(1);
opts = optiset('solver','cbc');
 Opt = opti('f',f,'ineq',A,b,'bounds',lb,ub,'xtype',xtype,'options',opts);
   [x,fval,exitflag,info] = solve(Opt);
% serverPos=find(x(1:probSize(1))>0.9);
% flow1=x(probSize(1)+1:probSize(1)+probSize(2));
% flow2=x(probSize(1)+probSize(2)+1:probSize(1)+probSize(2)+probSize(2));
% final=[edges,flow1,flow2];
% 
%  for i=1:96
%      if(final(i,6)~=0 || final(i,5)~=0)
%          display([final(i,1),final(i,2),final(i,5),final(i,6)]);
%      end
%  end
% flowMat=zeros(probSize(1),probSize(1));
% for i=1:50
%     for j=1:50
%         if(mat(i,j)>0)
%             flowMat(i,j)=final(mat(i,j),5);
%             flowMat(j,i)=final(mat(i,j),6);
%         end
%     end
% end