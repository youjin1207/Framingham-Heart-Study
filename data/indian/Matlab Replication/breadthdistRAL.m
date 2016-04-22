function  [R,D] = breadthdistRAL(CIJ,dummies)

% input:  
%           CIJ = connection/adjacency matrix
%           dummies = Dummy vector of length N for which you wish to
%                       compute the distance and reachability matrices.
% outputs: 
%           R   = reachability matrix
%           D   = distance matrix

% This function is potentially less memory-hungry than 'reachdist.m',
% particularly if the characteristic path length is rather large.
%
% Olaf Sporns, Indiana University, 2002/2007/2008

N = size(CIJ,1);

% D = zeros(N,sum(dummies));
D = zeros(N); % I'm not sure why I used the command above. It doesn't actually do what I think it was intended to do.
for i=find(dummies)',
   D(:,i) = breadth(CIJ,i);
end;

% replace zeros with 'Inf's
D(isinf(D)) = 999999;

% construct R
R = double(D~=999999);

