function [infectedbefore, dynamicInfection, contagious] = diffusion_InformationModel(parms,Z,Betas,X,leaders,j,T, EmpRate)

qN = parms(1); % Probability non-taker transmits information.
qP = parms(2); % Probability that a just-informed-taker transmits information.
N = size(X,1); % Number of individuals.

infected = false(N,1); % Nobody has been infected yet.
infectedbefore = false(N,1); % Nobody has been infected yet.
contagiousbefore = false(N,1); % People who were contagious before
contagious = logical(leaders); % Newly informed/contagious.
dynamicInfection = []; % Will be a vector that tracks the infection rate for the number of periods it takes place


x = rand(N,T);
t = 1;
for t = 1:T
    qNt = qN; % We may want to make qNt a dynamic function of qN
    qPt = qP; % We may want to make qNt a dynamic function of qP

    % Step 1: Take-up decision based on Newly informed, p.
    LOGITprob = 1./(1+exp(-([ones(N,1) Z]*Betas)));
    infected = ((~contagiousbefore & contagious & x(:,t) < LOGITprob)  | infected);
    s1 = sum(infected);
    s2 = sum(infectedbefore);
    infectedbefore = (infectedbefore | infected);
    contagiousbefore = (contagious | contagiousbefore);
    C = sum(contagious);

    % Step 2: Information flows
    transmitPROB = (contagious & infected)*qPt + (contagious & ~infected)*qNt; % Probability of transmission: Pr(.|infected)==qN, Pr(.|~infected)==qP. Individual (node) specific.
    contagionlikelihood = X(contagious,:).*(transmitPROB(contagious)*ones(1,N));

    % Step 3:
    contagious = ( ((contagionlikelihood > rand(C,N))'*ones(C,1) > 0) | contagiousbefore);
%     [t j s1 s2 sum(infectedbefore)/N EmpRate];
    dynamicInfection = [dynamicInfection; sum(infectedbefore)/N];
    t = t+1;
    [qN qP t];
end

