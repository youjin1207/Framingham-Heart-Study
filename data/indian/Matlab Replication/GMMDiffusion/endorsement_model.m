function [infectedbeforeE, infectedbeforeD, infectedbeforeN, dynamicInfectionE, dynamicInfectionD, dynamicInfectionN, contagiousbeforeE, contagiousbeforeD, contagiousbeforeN] = endorsement_model(parms,Z,Betas,X,leaders,OmegaE, OmegaD, OmegaN, j, T, EmpRate)
% Diffusion Model with three parameters.

% 3 Parameters: qN, qP, lambda.
qN = parms(1); % Probability non-taker transmits information.
qP = parms(2); % Probability that a just-informed-taker transmits information.
lambda = parms(3); % Coefficient in logit governing the probability of being a taker.

% Programmed by Arun Chandrasekhar, Nov 2010
% Originally "diffusionmodelNewneighborCovXval.m"
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(X,1); % Number of individuals.

% Eigenvector centrality
infectedE = false(N,1); % Nobody has been infected yet.
infectedbeforeE = false(N,1); % Nobody has been infected yet.
contagiousbeforeE = false(N,1); % Nobody has been contagious yet.
contagiousE = logical(leaders); % Newly informed/contagious.
transmissionHistE = false(N,N); % Who has transmitted to whom?
dynamicInfectionE = []; % Will be a vector that tracks the infection rate for the number of periods it takes place

% Degree
infectedD = false(N,1); % Nobody has been infected yet.
infectedbeforeD = false(N,1); % Nobody has been infected yet.
contagiousbeforeD = false(N,1); % Nobody has been contagious yet.
contagiousD = logical(leaders); % Newly informed/contagious.
transmissionHistD = false(N,N); % Who has transmitted to whom?
dynamicInfectionD = []; % Will be a vector that tracks the infection rate for the number of periods it takes place

% Naive
infectedN = false(N,1); % Nobody has been infected yet.
infectedbeforeN = false(N,1); % Nobody has been infected yet.
contagiousbeforeN = false(N,1); % Nobody has been contagious yet.
contagiousN = logical(leaders); % Newly informed/contagious.
transmissionHistN = false(N,N); % Who has transmitted to whom?
dynamicInfectionN = []; % Will be a vector that tracks the infection rate for the number of periods it takes place


x = rand(N,T);

for t = 1:T
    % Step Pre-1: Take-up decision
    %construct regressor
    zE = infectedbeforeE*ones(1,N);
    regressorE = diag(OmegaE*(transmissionHistE.*zE))./diag(OmegaE*(transmissionHistE));
    regressorE(isnan(regressorE))=0;
    regressorE(isinf(regressorE))=0;
    LOGITprobE = 1./(1+exp(-([ones(N,1) Z]*Betas + regressorE*lambda)));
    
    infectedE = ((~contagiousbeforeE & contagiousE & x(:,t) < LOGITprobE) | infectedE);
    s1 = sum(infectedE);
    s2 = sum(infectedbeforeE);
    infectedbeforeE = (infectedbeforeE | infectedE);
    contagiousbeforeE = (contagiousE | contagiousbeforeE);
    CE = sum(contagiousE);
    
    % Step 2: Information flows
    transmitPROBE = (contagiousE & infectedE)*qP + (contagiousE & ~infectedE)*qN; % Probability of transmission: Pr(.|infected)==qN, Pr(.|~infected)==qP. Individual (node) specific.
    contagionlikelihoodE = X.*(transmitPROBE*ones(1,N)); % A full NxN matrix of transmission probabilities
    
    % Step 3:
    t0E = rand(N,N);
    t1E = (contagionlikelihoodE > t0E); % a full transmission matrix
    
    % zero out stuff because one can only transmit if contagious, and one
    % cannot be transmitted to unless they were not contagious before
    t1E(contagiousE,contagiousbeforeE)=0; % This does not need to be done, but it makes no difference.
    t1E(~contagiousE,~contagiousbeforeE)=0;
    t1E(~contagiousE,contagiousbeforeE)=0;
    
    transmissionHistE = (transmissionHistE | t1E);
    
    t2E = t1E(contagiousE,~contagiousbeforeE); % which contagious folk transmit to previously uncontagious
    contagiousE(~contagiousbeforeE) = ( (t2E'*ones(CE,1) > 0));
    
    dynamicInfectionE = [dynamicInfectionE; sum(infectedbeforeE)/N];
end




for t = 1:T
    % Step Pre-1: Take-up decision
    %construct regressor
    zD = infectedbeforeD*ones(1,N);
    regressorD = diag(OmegaD*(transmissionHistD.*zD))./diag(OmegaD*(transmissionHistD));
    regressorD(isnan(regressorD))=0;
    regressorD(isinf(regressorD))=0;
    LOGITprobD = 1./(1+exp(-([ones(N,1) Z]*Betas + regressorD*lambda)));
    
    infectedD = ((~contagiousbeforeD & contagiousD & x(:,t) < LOGITprobD) | infectedD);
    s1 = sum(infectedD);
    s2 = sum(infectedbeforeD);
    infectedbeforeD = (infectedbeforeD | infectedD);
    contagiousbeforeD = (contagiousD | contagiousbeforeD);
    CD = sum(contagiousD);
    
    % Step 2: Information flows
    transmitPROBD = (contagiousD & infectedD)*qP + (contagiousD & ~infectedD)*qN; % Probability of transmission: Pr(.|infected)==qN, Pr(.|~infected)==qP. Individual (node) specific.
    contagionlikelihoodD = X.*(transmitPROBD*ones(1,N)); % A full NxN matrix of transmission probabilities
    
    % Step 3:
    t0D = rand(N,N);
    t1D = (contagionlikelihoodD > t0D); % a full transmission matrix
    
    % zero out stuff because one can only transmit if contagious, and one
    % cannot be transmitted to unless they were not contagious before
    t1D(contagiousD,contagiousbeforeD)=0;
    t1D(~contagiousD,~contagiousbeforeD)=0;
    t1D(~contagiousD,contagiousbeforeD)=0;
    
    transmissionHistD = (transmissionHistD | t1D);
    
    t2D = t1D(contagiousD,~contagiousbeforeD); % which contagious folk transmit to previously uncontagious
    %     contagiousD = false(N,1);
    contagiousD(~contagiousbeforeD) = ( (t2D'*ones(CD,1) > 0));
    
    dynamicInfectionD = [dynamicInfectionD; sum(infectedbeforeD)/N];
end




for t=1:T
    % Step Pre-1: Take-up decision
    %construct regressor
    zN = infectedbeforeN*ones(1,N);
    regressorN = diag(OmegaN*(transmissionHistN.*zN))./diag(OmegaN*(transmissionHistN));
    regressorN(isnan(regressorN))=0;
    regressorN(isinf(regressorN))=0;
    LOGITprobN = 1./(1+exp(-([ones(N,1) Z]*Betas + regressorN*lambda)));
    
    infectedN = ((~contagiousbeforeN & contagiousN & x(:,t) < LOGITprobN) | infectedN);
    s1 = sum(infectedN);
    s2 = sum(infectedbeforeN);
    infectedbeforeN = (infectedbeforeN | infectedN);
    contagiousbeforeN = (contagiousN | contagiousbeforeN);
    CN = sum(contagiousN);
    
    % Step 2: Information flows
    transmitPROBN = (contagiousN & infectedN)*qP + (contagiousN & ~infectedN)*qN; % Probability of transmission: Pr(.|infected)==qN, Pr(.|~infected)==qP. Individual (node) specific.
    contagionlikelihoodN = X.*(transmitPROBN*ones(1,N)); % A full NxN matrix of transmission probabilities
    
    % Step 3:
    t0N = rand(N,N);
    t1N = (contagionlikelihoodN > t0N); % a full transmission matrix
    
    % zero out stuff because one can only transmit if contagious, and one
    % cannot be transmitted to unless they were not contagious before
    t1N(contagiousN,contagiousbeforeN)=0;
    t1N(~contagiousN,~contagiousbeforeN)=0;
    t1N(~contagiousN,contagiousbeforeN)=0;
    
    transmissionHistN = (transmissionHistN | t1N);
    
    t2N = t1N(contagiousN,~contagiousbeforeN); % which contagious folk transmit to previously uncontagious
    contagiousN(~contagiousbeforeN) = ( (t2N'*ones(CN,1) > 0));
    
    dynamicInfectionN = [dynamicInfectionN; sum(infectedbeforeN)/N];
end
