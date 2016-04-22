%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIVERGENCE_Model
% This computes the deviation of the empirical moments from the simulated
% ones for:

% Model 1 = q
% Model 3 = qN, qP
% It relies on diffusion_model as the transmission process and moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D TimeSim] = divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, T,EmpRate, version)

%% Parameters
G = length(X);


%% Computation of the vector of divergences across all the moments
EmpiricalMoments = zeros(G,m);
MeanSimulatedMoments = zeros(G,m);
D = zeros(G,m);
TimeSim = zeros(G,S);

for g=1:G
    % Compute moments - G x m object
    EmpiricalMoments(g,:) = moments(X{g},leaders{g},TakeUp{g},Sec{g},g,version);
    
    % Compute simulated moments
    SimulatedMoments = zeros(S,m);
    for s=1:S
        infectedSIM = diffusion_model(theta, Z{g}, Betas, X{g},leaders{g}, g, T(g), EmpRate(g));
        SimulatedMoments(s,:) = moments(X{g},leaders{g},infectedSIM,Sec{g},g,version);
    end    
    
    % Compute the mean simulated moment - a G x m object
    MeanSimulatedMoments(g,:) = mean(SimulatedMoments,1);   
    D(g,:) = MeanSimulatedMoments(g,:) - EmpiricalMoments(g,:);
    ['Done with ' num2str(g/G*100) '% of the graphs for THIS parameter value']
    theta
end

