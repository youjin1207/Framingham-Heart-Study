%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Social Networks and Microfinance: GMM Bootstrap
%
%
% Programmed by Arun Chandrasekhar May 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DE, DD, DN, TSimsE, TSimsD, TSimsN] = divergence_endorsement_model(X, Z, Betas, leaders, TakeUp, OmegaE, OmegaD, OmegaN, Sec, theta, m, S, T, EmpRate, version)


%% Parameters
G = length(X); % Number of graphs

%% Computation of the vector of divergences across all the moments
EmpiricalMoments = zeros(G,m);
MeanSimulatedMomentsE = zeros(G,m);
MeanSimulatedMomentsD = zeros(G,m);
MeanSimulatedMomentsN = zeros(G,m);
DE = zeros(G,m);
DD = zeros(G,m);
DN = zeros(G,m);
% TSimsE = zeros(G,S);
% TSimsD = zeros(G,S);
% TSimsN = zeros(G,S);

for g=1:G
    % Compute moments - G x m object
    EmpiricalMoments(g,:) = moments(X{g},leaders{g},TakeUp{g},Sec{g},g,version);
    % Compute simulated moments
    SimulatedMomentsE = zeros(S,m);
    SimulatedMomentsD = zeros(S,m);
    SimulatedMomentsN = zeros(S,m);
    for s=1:S
        oneSimTime = tic;
        [InfectedSimsE, InfectedSimsD, InfectedSimsN, TSimsE, TSimsD, TSimsN] = endorsement_model(theta, Z{g}, Betas, X{g},leaders{g}, OmegaE{g}, OmegaD{g}, OmegaN{g}, g, T(g), EmpRate(g));
        SimulatedMomentsE(s,:) = moments(X{g},leaders{g},InfectedSimsE,Sec{g},g,version);
        SimulatedMomentsD(s,:) = moments(X{g},leaders{g},InfectedSimsD,Sec{g},g,version);
        SimulatedMomentsN(s,:) = moments(X{g},leaders{g},InfectedSimsN,Sec{g},g,version);
        toc(oneSimTime)
        [s g theta]
        ['Done with ' num2str(g/G*100) '% of the graphs for THIS parameter value']
    end

    % Compute the mean simulated moment - a G x m object
    MeanSimulatedMomentsE(g,:) = mean(SimulatedMomentsE,1);
    MeanSimulatedMomentsD(g,:) = mean(SimulatedMomentsD,1);
    MeanSimulatedMomentsN(g,:) = mean(SimulatedMomentsN,1);
    
    DE(g,:) = MeanSimulatedMomentsE(g,:) - EmpiricalMoments(g,:);
    DD(g,:) = MeanSimulatedMomentsD(g,:) - EmpiricalMoments(g,:);
    DN(g,:) = MeanSimulatedMomentsN(g,:) - EmpiricalMoments(g,:);
end

