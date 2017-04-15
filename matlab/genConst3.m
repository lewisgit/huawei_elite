function [f,intcon,A,b,Aeq,beq,lb,ub]=genConst3(probSize,servCost,edges,conReq,mat,x)
var1=probSize(1);
var2=probSize(2);
conNum=size(conReq,1);
varNum=var1+var2*2;
%% compute f;
f=ones(1,var1)*servCost;

%additional edges
% additions=zeros(1,conNum);
% additions=[edges(:,4)',additions];
f=[f,edges(:,4)',edges(:,4)',zeros(1,2*var1+1)];
% f=-f;
%% compute intcon;
intcon=1:var1;
%% compute Aeq;
Aeq=[];
beq=[];
%% compute A;
subMat1=eye(var1)*-100000;
subMat2=zeros(var1,var2);
for i=1:var1
    for j=1:var1
        if(mat(i,j)>0)
            subMat2(i,mat(i,j))=1;
        end
        if(mat(j,i)>0)
            subMat2(i,mat(j,i))=-1;
        end
        %tmp=find(conReq(:,2)==(i-1));

    end
end
subMat2=[subMat2,-subMat2];

% for i=1:var1
%     for j=1:var1
%         tmp=find(conReq(:,2)==(i-1));
%         if(~isempty(tmp))
%             subMat2(i,var2+tmp)=conReq(tmp,3);
%         end
%     end
% end


A=[subMat1,subMat2;
    subMat1,-subMat2];
% subMat2=zeros(conNum,var2);
% for i=1:conNum
%     for j=1:var1
%         if(mat(conReq(i,2)+1,j)>0)
%             subMat2(i,mat(conReq(i,2)+1,j))=-1;
%         end
%     end
% end
%subMat2=[zeros(conNum,var1),subMat2,-subMat2];
%A=[A;subMat2];
A=[A;ones(1,var1),zeros(1,2*var2)];
% A=[A;eye(varNum)];
A=[A,eye(2*var1+1)];
% A=[A,[zeros(2*var1,1);1]];
% for i=1:conNum
%     A(conReq(i,2)+1,:)=zeros(1,varNum);
%     A(conReq(i,2)+1+var1,:)=zeros(1,varNum);
% end

b=zeros(var1*2+1,1);
%b=[b;-conReq(:,3)];
% b=[b;conNum];

for i=1:size(conReq,1)
    b(conReq(i,2)+1)=-conReq(i,3);
    b(conReq(i,2)+1+var1)=conReq(i,3);
%     tmp=find(conReq(:,2)==(i-1));
%     if(~isempty(tmp))
%         b(i)=-conReq(tmp,3);
%         b(var1+i)=conReq(tmp,3);
%     end
end
b(var1*2+1)=conNum;
% b(var1*2+2:var1*3+1)=1;
% b(var1*3+2:2*var1+varNum+1)=[edges(:,3);edges(:,3)];

lb=zeros(varNum+2*var1+1,1);

% 
ub=[ones(var1,1);edges(:,3);edges(:,3);inf*ones(2*var1+1,1)];
% lb(1:50)=x(1:50);
% ub(1:50)=x(1:50);
% for i=var1+var2+1:var1+var2+conNum;
%     lb(i)=1;
%     ub(i)=1;
%     lb(i+var2+conNum)=0;
%     ub(i+var2+conNum)=0;
% end