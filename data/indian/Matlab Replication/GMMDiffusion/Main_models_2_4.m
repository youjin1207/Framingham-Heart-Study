%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main_models_2_4

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

tic;
% Adjusting so this runs off dropbox
cd ..
clear all

location = pwd;
addpath(genpath(location));

%% Parameters
% 0. Model type
modelType = 4; % modelType 2 if q, modelType 4 if qN and qP, qN != qP

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
        T = ceil(TMonths./3) + 1  %Quarters have 3 months in them
    case 'trimesters'
        T = ceil(TMonths./4) + 1 %Trimesters have 4 months in them
end
assert(G == numel(T))


% 4. Select parameter grid
lambda = [(-1 : .1 : -0.3), (-0.25 : .05 : 0.3), (0.4 : .1 : 1)];

if modelType == 2,
    qN = [(0 : 0.05 : 0.5), (0.6 : 0.1: 1)];
elseif modelType == 4,
    qN = [(0 : 0.05 : 0.5), (0.6 : 0.1: 1)];
    qP = [(0 : 0.05 : 0.5), (0.6 : 0.1: 1)];
end

% 5. Load data
relative = 1; 

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
Omega_E = cell(G, 1);
Omega_D = cell(G, 1);
Omega_N = cell(G, 1);


% Load the household connection adjacency matrix.
X = load(['India Networks/adjacencymatrix.mat']);
X = X.X;

%% Construct data
counter = 0;
for vilnum = vills
    counter = counter + 1;
   
    N = length(X{counter});
    
    if relative==0
        Omega_E{counter} = dlmread(['./India Networks/Omega_abs' num2str(vilnum) '.csv']);
        Omega_D{counter} = dlmread(['./India Networks/DOmega_abs' num2str(vilnum) '.csv']);
        Omega_N{counter} = dlmread(['./India Networks/NOmega_abs' num2str(vilnum) '.csv']);
    elseif relative==1,
        Omega_E{counter} = dlmread(['./India Networks/Omega_rel' num2str(vilnum) '.csv']);
        Omega_D{counter} = dlmread(['./India Networks/DOmega_rel' num2str(vilnum) '.csv']);
        Omega_N{counter} = dlmread(['./India Networks/NOmega_rel' num2str(vilnum) '.csv']);
    end
    
    
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
    Omega_E{counter} = Omega_E{counter}(logical(inGiant{counter}),logical(inGiant{counter}));
    Omega_D{counter} = Omega_N{counter}(logical(inGiant{counter}),logical(inGiant{counter}));
    Omega_N{counter} = Omega_N{counter}(logical(inGiant{counter}),logical(inGiant{counter}));
    
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


% 6. Logistic fit to get coefficients for covariates and the constant
[Betas, dev, stats] = glmfit(Covars,Outcome,'binomial','link','logit');
[Betas'; stats.se'; stats.p'];

tic;

% 7. RUNNING THE MODEL
%% Obtain the divergence matrix
if modelType==2, % case where qN = qP
    %% Obtain the divergence matrix
    DE = cell(length(qN),length(lambda));
    DD = cell(length(qN),length(lambda));
    DN = cell(length(qN),length(lambda));
    TimeSimE = cell(length(qN),length(lambda));
    TimeSimD = cell(length(qN),length(lambda));
    TimeSimN = cell(length(qN),length(lambda));
    iterations = 0;
    totalCount = length(qN)*length(lambda);
    for i=1:length(qN)
        for j=1:length(lambda)
            
            oneGridPtTime = tic;
            theta = [qN(i), qN(i), lambda(j)]
            [DE{i,j}, DD{i,j}, DN{i,j}, TimeSimE{i,j}, TimeSimD{i,j}, TimeSimN{i,j}] = divergence_endorsement_model(X, Z, Betas, leaders, TakeUp, Omega_E, Omega_D, Omega_N, Sec, theta, m, S, T, EmpRate, version);
            
            iterations = iterations+1;
            ['Done with ' num2str(iterations/totalCount*100) '% of the D(G,m) computations.']
            oneGridPtTimeelasped = toc(oneGridPtTime);
            
        end
    end
    
elseif modelType==4, % case where qN and qP are independent
    %% Obtain the divergence matrix
    DE = cell(length(qN), length(qP), length(lambda));
    DD = cell(length(qN),length(qP), length(lambda));
    DN = cell(length(qN),length(qP), length(lambda));
    TimeSimE = cell(length(qN),length(qP), length(lambda));
    TimeSimD = cell(length(qN),length(qP),length(lambda));
    TimeSimN = cell(length(qN),length(qP),length(lambda));
    iterations = 0;
    totalCount = length(qN)*length(qP)*length(lambda);
    for i=1:length(qN)
        for j=1:length(qP)
            for k = 1:length(lambda)
                oneGridPtTime = tic;
                theta = [qN(i), qP(j), lambda(k)]
                [DE{i,j,k}, DD{i,j,k}, DN{i,j,k}, TimeSimE{i,j,k}, TimeSimD{i,j,k}, TimeSimN{i,j,k}] = divergence_endorsement_model(X, Z, Betas, leaders, TakeUp, Omega_E, Omega_D, Omega_N, Sec, theta, m, S, T, EmpRate, version);
                iterations = iterations+1;
                ['Done with ' num2str(iterations/totalCount*100) '% of the D(G,m) computations.']
                oneGridPtTimeelasped = toc(oneGridPtTime);
            end
        end
        
    end
end
toc;

% Process Data - construct data names
outputName = ['data_model_' num2str(modelType) '_mom_' num2str(version) '']
save([outputName, ' ', timeVector, ' ','.mat'])

% 8. RUNNING THE AGGREGATOR
% Bootstrap?
bootstrap=0; % yes or no
if bootstrap==0,
    B = 1;
elseif bootstrap==1
    B = 1000;
end

%% Set up data
if model ==2,
    q = qN;
    % Set up matrices
    DivDD = zeros(length(q), length(lambda), G, m);
    DivDE = zeros(length(q), length(lambda), G, m);
    DivDN = zeros(length(q), length(lambda), G, m);
elseif model == 4,
    % Set up matrices
    DivDD = zeros(length(qN), length(qP), length(lambda), G, m);
    DivDE = zeros(length(qN), length(qP), length(lambda), G, m);
    DivDN = zeros(length(qN), length(qP), length(lambda), G, m);
end

%% Put data into matrix
if model == 2,
    for i = 1:length(q)
        for k = 1:length(lambda)
            DDTotal(i, k, :, :) = DD{i, k};
            DETotal(i, k, :, :) = DE{i, k};
            DNTotal(i, k, :, :) = DN{i, k};
        end
    end
elseif model == 4,
    for i = 1:length(qN)
        for j = 1:length(qP)
            for k = 1:length(lambda)
                DDTotal(i,j,k, :, :) = DD{i, j,k};
                DETotal(i, j, k, :, :) = DE{i, j, k};
                DNTotal(i, j,k, :, :) =DN{i,j,k};
            end
        end
    end
    
end

%Two Step Optimal Weights
twoStepOptimal = 0;

if twoStepOptimal == 1,
    % [qN_E, qP_E, lambda_E] = ;
    theta = [qN_E, qP_E, lambda_E];
    [DE, DD, DN, TSimsE, TsimsD, TsimsN]= divergence_endorsement_model(X, Z, Betas, leaders, TakeUp,Omega_E, Omega_D, Omega_N, Sec, theta, m, S, T, EmpRate, version);
    AE = (DE'*DE)/43;
    WE = AE^(-1);
    
    % [qN_D, qP_D, lambda_D] = ;
    theta = [qN_D, qP_D, lambda_D];
    [DE, DD, DN, TSimsE, TsimsD, TsimsN]= divergence_endorsement_model(X, Z, Betas, leaders, TakeUp,Omega_E, Omega_D, Omega_N, Sec, theta, m, S, T, EmpRate, version);
    AD = (DD'*DD)/43;
    WD = AD^(-1);
    
    % [qN_N, qP_N, lambda_N] = ;
    theta = [qN_N, qP_N, lambda_N];
    [DE, DD, DN, TSimsE, TsimsD, TsimsN]= divergence_endorsement_model(X, Z, Betas, leaders, TakeUp,Omega_E, Omega_D, Omega_N, Sec, theta, m, S, T, EmpRate, version);
    AN = (DN'*DN)/43;
    WN = AN^(-1);    
else
    WE = eye(m);
    WD = eye(m);
    WN = eye(m);
end


% Pre-allocation
QEndE = zeros(B,1);
QEndD = zeros(B,1);
QEndN = zeros(B,1);
TestEndEEndD = zeros(B,1);
TestEndDEndN = zeros(B,1);
TestEndEEndN = zeros(B,1);
importantparmsE = [];
importantparmsD = [];
importantparmsN = [];
valE = [];
valD = [];
valN = [];


% Aggregate
for b=1:B
    % Generate weights b
    if bootstrap==1,
        wt(b,:) = exprnd(1,G,1);
        wt(b,:) = wt(b,:)/mean(wt(b,:));
    elseif bootstrap==0
        wt(b,:) = 1/G*ones(1,G);
    end
    
    
    % For each model, generate the criterion function value for this
    % bootstrap run
    
    % ENDORSEMENT MODEL
    for i=1:length(qN)
        for j=1:length(qP)
            for k=1:length(lambda)
                % Compute the moment function
                momFuncE(i,j, k,b,:) = (wt(b,:)*squeeze(DETotal(i,j, k,:,:)))'/G;
                momFuncD(i,j, k,b,:) = (wt(b,:)*squeeze(DDTotal(i,j, k,:,:)))'/G;
                momFuncN(i,j, k,b,:) = (wt(b,:)*squeeze(DNTotal(i,j, k,:,:)))'/G;
                
                
                % Criterion function
                QEndorseE(i,j, k,b) = (squeeze(momFuncE(i,j, k,b,:)))'*WE*squeeze(momFuncE(i,j, k,b,:));
                QEndorseD(i,j, k,b) = (squeeze(momFuncD(i,j, k,b,:)))'*WD*squeeze(momFuncD(i,j, k,b,:));
                QEndorseN(i,j, k,b) = (squeeze(momFuncN(i,j, k,b,:)))'*WN*squeeze(momFuncN(i,j, k,b,:));
                
            end
        end
        
    end
    
    TempEndE = QEndorseE(:,:, :,b);
    TempEndD = QEndorseD(:,:, :,b);
    TempEndN = QEndorseN(:,:, :,b);
    
    [minAEndE,indEndE] = min(TempEndE(:));
    [x1EndE,x2EndE, x3EndE] = ind2sub(size(TempEndE),indEndE);
    [minAEndD,indEndD] = min(TempEndD(:));
    [x1EndD,x2EndD, x3EndD] = ind2sub(size(TempEndD),indEndD);
    [minAEndN,indEndN] = min(TempEndN(:));
    [x1EndN, x2EndN, x3EndN] = ind2sub(size(TempEndN),indEndN);
    
    importantparmsE = [importantparmsE; x1EndE,x2EndE, x3EndE];
    importantparmsD = [importantparmsD; x1EndD,x2EndD, x3EndD];
    importantparmsN = [importantparmsN; x1EndN,x2EndN, x3EndN];
    valE = [valE; qN(x1EndE), qP(x2EndE), lambda(x3EndE)];
    valD = [valD; qN(x1EndD), qP(x2EndD), lambda(x3EndD)];
    valN = [valN; qN(x1EndN), qP(x2EndN), lambda(x3EndN)];
    
    % Need to map back
    QEndE(b) = QEndorseE(x1EndE,x2EndE, x3EndE, b);
    QEndD(b) = QEndorseD(x1EndD,x2EndD, x3EndD, b);
    QEndN(b) = QEndorseN(x1EndN,x2EndN, x3EndN, b);
    
    % Tests
    TestEndEEndD(b) = sqrt(G)*(QEndE(b)-QEndD(b));
    TestEndDEndN(b) = sqrt(G)*(QEndD(b)-QEndN(b));
    TestEndEEndN(b) = sqrt(G)*(QEndE(b)-QEndN(b));
    
    ['Done with ' num2str(100*b/B) '% of the bootstraps']
    
end

%% The test
CIendEEndD = [prctile(TestEndEEndD,5), prctile(TestEndEEndD,50), prctile(TestEndEEndD,95)]
CIendDEndN = [prctile(TestEndDEndN,5), prctile(TestEndDEndN,50), prctile(TestEndDEndN,95)]
CIendEEndN = [prctile(TestEndEEndN,5), prctile(TestEndEEndN,50), prctile(TestEndEEndN,95)]

[mean(valE); std(valE); median(valE)]
[mean(valD); std(valD); median(valD)]
[mean(valN); std(valN); median(valN)]

diffqNqPE = valE(:,1) - valE(:,2);
diffqNqPD = valD(:,1) - valD(:,2);
diffqNqPN = valN(:,1) - valN(:,2);

[mean(diffqNqPE); std(diffqNqPE)]
[mean(diffqNqPD); std(diffqNqPD)]
[mean(diffqNqPN); std(diffqNqPN)]

