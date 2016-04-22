%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main_models_1_3

% Outline

% 0. Model type
% 1. Select villages to consider
% 2. Select moments
% 3. Select time vector and number of repetitions per trial
% 4. Select parameter grid
% 5. Load data
% 6. Logistic fit to get coefficients for covariates and the constant
% 7. RUNNING THE MODEL
% 8. RUNNING THE AGGREGATOR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
tic;
% Adjusting so this runs off dropbox
cd ..
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


% 4. Select parameter grid
if modelType == 1,
    qN = [(0:0.001:0.01), (0.05:0.05:1)];
elseif modelType == 3,
    qN = [(0:0.001:0.01), (0.05:0.05:1)];
    qP = [(0:0.005:0.1), (0.15:0.05:1)];
end


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


%7. RUNNING THE MODEL
if modelType==1, % case where qN = qP
    D = cell(length(qN));
    TimeSim = cell(length(qN));
    iterations = 0;
    totalCount = length(qN);
    
    for i=1:length(qN)
        theta = [qN(i), qN(i)];
        oneGridPtTime = tic;
        [D{i} TimeSim{i}] = divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, T, EmpRate, version);
        iterations = iterations+1;
        ['Done with ' num2str(iterations/totalCount*100) '% of the D(G,m) computations.']
        oneGridPtTimeElasped=toc(oneGridPtTime);
    end
    
    
elseif modelType==3, % case where qN and qP are independent
    D = cell(length(qN),length(qP));
    TimeSim = cell(length(qN),length(qP));
    iterations = 0;
    totalCount = length(qN)*length(qP);
    
    for i=1:length(qN)
        for j=1:length(qP)
            theta = [qN(i), qP(j)];
            [D{i,j} TimeSim{i,j}] = divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, T, EmpRate, version);
            iterations = iterations+1;
            ['Done with ' num2str(iterations/totalCount*100) '% of the D(G,m) computations.']
            toc;
        end
    end
end


% Process Data - construct data names
outputName = ['data_model_' num2str(modelType) '_mom_' num2str(version) '']
save([outputName, ' ', timeVector, ' ','.mat'])


%8. RUNNING THE AGGREGATOR
% Bootstrap?
bootstrap=0;
if bootstrap==0,
    B = 1;
elseif bootstrap==1
    B = 1000;
end
Q = zeros(B,1);

%Two Step Optimal Weights
twoStepOptimal = 0;
if twoStepOptimal == 1,
    qN_info = 0.09;
    qP_info = 0.45;
    theta = [qN_info, qP_info];
    D = divergence_model(X, Z, Betas, leaders, TakeUp, Sec, theta, m, S, T, EmpRate, version);
    A = (D'*D)/43;
    W = A^(-1);
else
    W = eye(m);
end


%Put data into matrix
Dnew = zeros(length(qN), length(qP), G, m);
for i = 1:length(qN)
    for j = 1:length(qP)
        Dnew(i, j, :, :) =D{i,j};
    end
end

% Pre-allocation
importantparms = [];
val = [];

% Aggregate
for b=1:B
    
    % Generate weights b  
    if bootstrap==1,
        wt(b,:) = exprnd(1,G,1);
        wt(b,:) = wt(b,:)/mean(wt(b,:));
    elseif bootstrap==0
        wt(b,:) = 1/G*ones(1,G);
    end
    
    %% For each model, generate the criterion function value for this
    %% bootstrap run
   
    % Info model
    for i=1:length(qN)
        for j=1:length(qP)
            % Compute the moment function
            momFunc(i,j,b,:) = (wt(b,:)*Dnew{i,j})'/G;        
            % Criterion function
            Qa(i,j,b) = (squeeze(momFunc(i,j,b,:)))'*W*squeeze(momFunc(i,j,b,:));
        end
    end
    
    Temp = Qa(:,:,b);
    
    [minA,ind] = min(Temp(:));
    [x1,x2] = ind2sub(size(Temp),ind);
    importantparms = [importantparms; x1,x2];
    val = [val; qN(x1), qP(x2)];
    
    % Need to map back
    Q(b) = Qa(x1, x2, b);
    
    
    ['Done with ' num2str(100*b/B) '% of the bootstraps']
    
end


%% The test
[mean(val); std(val)]
diffqNqP = val(:,1) - val(:,2);
[mean(diffqNqP); std(diffqNqP)]

toc;







