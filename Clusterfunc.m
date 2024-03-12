function [order,clusts]= Clusterfunc(CorrData,n)
dim=size(CorrData);
D=.5*ones(dim(1),dim(2));
D=D-.5*(CorrData);
clear CorrData
ZZ=linkage(squareform(D),'average');
order=optimalleaforder(ZZ,squareform(D),'Criteria','adjacent');
clusts=cluster(ZZ,'maxclust',n);
end