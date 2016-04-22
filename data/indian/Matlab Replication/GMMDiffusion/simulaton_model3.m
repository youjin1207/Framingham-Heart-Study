tic;
% Adjusting so this runs off dropbox
cd ..
clear all

location = pwd;
addpath(genpath(location));

%% Parameters
% 0. Model type
modelType = 3; % modelType 1 if q, modelType 3 if qN and qP, qN != qP

% 1. Select villages to consider
vills = [1:4,6,9, 12, 15, 19:21, 23:25, 29, 31:33, 36, 39, 42, 43, 45:48, 50:52, 55, 57, 59:60, 62, 64:65, 67:68, 70:73, 75];
G = length(vills); % Number of graphs

% 2. Select moments
version = 1;

switch version
    case 1
        m = 5;
    case 2
        m = 3;
    case 3
        m = 3;
    case 4
        m = 3;
    case 5
        m = 3;
end

% 3. Select time vector and number of repetitions per trial
S = 75; % Number of simulations
timeVector = 'trimesters'
TMonths = [31 35 15 35 13 2 32 5 31 35 31 29 19 22 25 25 23 23 24 25 26 24 17 16 17 13 19 20 20 19 22 14 12 15 10 19 18 18 19 19 19 17 17]; % Months

switch timeVector
    case 'months'
        T = TMonths + 1
    case 'quarters'
        T = ceil(TMonths./3) + 1 %Quarters have 3 months in them
    case 'trimesters'
        T = ceil(TMonths./4) + 1 %Trimesters have 4 months in them
end
assert(G == numel(T))

% 5. Load data
%% Pre-allocation
X = cell(G,1);
inGiant = cell(G,1);
leaders = cell(G,1);
TakeUp = cell(G,1);
TakingLeaders = cell(G,1);
ZLeaders = cell(G,1);
Covars = [];
Outcome = [];
EmpRate = zeros(G,1);
rdist = cell(G,1);
dist = cell(G,1);
hermits = cell(G,1);
W = cell(G,1);
Z = cell(G,1);
temp = cell(G, 1); 

% Load the household connection adjacency matrix.
X = load(['India Networks/adjacencymatrix.mat']);
X = X.X;

%% Construct data
counter = 0;
for vilnum = vills
    counter = counter + 1;
    
    N = length(X{counter});
    
    % Load the Leader data
    templeaders = load(['./India Networks/HHhasALeader' num2str(vilnum) '.csv']);
    leaders{counter} = templeaders(:,2);
    
    
    % Load the Take-Up data
    TakeUp{counter} = load(['./India Networks/MF' num2str(vilnum) '.csv']);
    EmpRate(counter) = mean(TakeUp{counter}(~leaders{counter}));
    
    % Load the giant component data
    inGiant{counter} = load(['./India Networks/inGiant' num2str(vilnum) '.csv']);
    
    % Generate hermits
    d = sum(X{counter},2);
    hermits{counter}=(d==0);
    
    % Load the Covariates data
    W{counter} = load(['./India Networks/hhcovariates' num2str(vilnum) '.csv']);
    
    % Which covariates to use - for instance we want to add a PCA
    Z{counter} = [W{counter}(:,1:6)]; % for instance, take the first covariate only
    
    % prune other stats
    leaders{counter} = leaders{counter}(logical(inGiant{counter}));
    TakeUp{counter} = TakeUp{counter}(logical(inGiant{counter}));
    Z{counter} = Z{counter}(logical(inGiant{counter}),:);
    
    TakingLeaders{counter} = TakeUp{counter}(logical(leaders{counter}));
    ZLeaders{counter} = Z{counter}(logical(leaders{counter}),:);
    Outcome = [Outcome; TakingLeaders{counter}];
    Covars = [Covars; ZLeaders{counter}];
    
    % Second neighbors
    Sec{counter} = (X{counter}^2>0);
    for i=1:length(X{counter})
        Sec{counter}(i,i) = 0;
    end
    Sec{counter}=(Sec{counter}-X{counter}>0);
end


%6. Logistic fit to get coefficients for covariates and the constant
[Betas, dev, stats] = glmfit(Covars,Outcome,'binomial','link','logit');
[Betas'; stats.se'; stats.p']



%% SIM MODEL
endogSwitch = 0; % 1 - turn on a simulation of the endogenous model instead
qN0 = 0.05;
qP0 = 0.35;

SimCount = 1000;
theta0 = [qN0, qP0]
dynamicOutcomes = [];
dynamicEndog = [];
infectionRate = zeros(G, SimCount);

if endogSwitch==1,
    for g=1:G
        percentileCell{g} = [1:1:T(g)]'/max(T(g));
    end
    
    for g=1:G
        for simcount = 1:SimCount
            [infectedSIM{g}, tEndog ,dynamicSIM{g,simcount}]=  diffusion_InformationModel_Endogenous(theta0, Z{g}, Betas, X{g},leaders{g}, g, T(g), EmpRate(g));
            
            % run a logistic regression of outcomes on time
            tEndog = tEndog-1;
            blogit(:,g,simcount) = glmfit([1:1:tEndog]',dynamicSIM{g,simcount},'binomial');
            
            % compute the values at they key percentiles
            keyTimes = tEndog*percentileCell{g};
            endogDynamicTakeTmp(:,simcount) = 1./(1+exp(-blogit(1,g,simcount) - blogit(2,g,simcount)*keyTimes));
            
        end
        endogDynamicTakeUp{g} = mean(endogDynamicTakeTmp,2);
        clear endogDynamicTakeTmp keyTimes blogit
        ['Done with ' num2str(g) ' of ' num2str(G) ' of the villages!']
    end
    
    for g=1:G
        dynamicEndog = [dynamicEndog; endogDynamicTakeUp{g}, ones(length(endogDynamicTakeUp{g}),1)*vills(g), [1:1:length(endogDynamicTakeUp{g})]'];
    end
    csvwrite('../dynamicEndogenous.csv',dynamicEndog);
    
    
else
    for simcount = 1:SimCount
        for g=1:G
            [infectedSIM{g} dynamicSIM{g,simcount}, contagiousSIM{g}]= diffusion_InformationModel(theta0, Z{g}, Betas, X{g},leaders{g}, g, T(g), EmpRate(g));
            infectionRate(g,simcount) = mean(infectedSIM{g}(logical((1-leaders{g}))));
            informedRate(g,simcount) = mean(contagiousSIM{g}(logical((1-leaders{g}))));
        end
    end
    infectionRateMean = mean(infectionRate,2)
    informedRateMean = mean(informedRate,2)
    
    
    for g=1:G
        for simcount = 1:SimCount
            dynamicData{g}(:,simcount) = dynamicSIM{g,simcount};
        end
        dynamicMean{g} = mean(dynamicData{g},2);
    end
    
    for g=1:G
        dynamicOutcomes = [dynamicOutcomes; dynamicMean{g}, ones(length(dynamicMean{g}),1)*vills(g), [1:1:length(dynamicMean{g})]'];
    end
    csvwrite('../dynamicOutcomesTrim_trimesters_1000sim.csv',dynamicOutcomes);
    
end
toc

