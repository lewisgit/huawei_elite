function [mat]=genMat(edges,n)
mat=zeros(n,n);
for i=1:size(edges,1);
    mat(edges(i,1)+1,edges(i,2)+1)=i;
    %mat(edges(i,2)+1,edges(i,1)+1)=i;
end
