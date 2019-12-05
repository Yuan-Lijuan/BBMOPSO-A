function pop = initialize(popsize, bounds,k)
%realnumber GA initialize
numVars=size(bounds,1);%The number of variables
range=(bounds(:,2)-bounds(:,1))';%The variable ranges
pop=zeros(popsize,numVars);
pop(:,1:numVars)=(ones(popsize,1)*bounds(:,1)')+(ones(popsize,1)*range).*(rand(popsize,numVars));%精华之处**!!
%pop(:,numVars+1)=rand(popsize,1);


