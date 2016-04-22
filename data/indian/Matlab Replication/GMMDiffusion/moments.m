function stats = moments(X,leaders,infected, Sec,j, version)

persistent netstats

N = size(X,1);

if isempty(netstats) | (length(netstats) < j),
    [R,D] = breadthdistRAL(X, leaders);
    
    minDistFromLeaders = min(D(:,logical(leaders))')';
    avgDistFromLeaders = mean(D(:,logical(leaders))')';
    
    if sum(infected.*leaders)>0
        minDistInfectedLeaders = min(D(:,logical(infected.*leaders))')';
    else
        minDistInfectedLeaders = 0;
    end
    if sum((1-infected).*leaders)>0
        minDistNonInfectedLeaders = min(D(:,logical((1-infected).*leaders))')';
    else
        minDistNonInfectedLeaders = 0;
    end
    
    netstats(j).minDistFromLeaders = minDistFromLeaders;
    netstats(j).avgDistFromLeaders = avgDistFromLeaders;
    netstats(j).minDistInfectedLeaders = minDistInfectedLeaders;
    netstats(j).minDistNonInfectedLeaders = minDistNonInfectedLeaders;
    
    [size(netstats(j).minDistFromLeaders) j];
    [size(netstats(j).minDistInfectedLeaders) j];
    [size(netstats(j).minDistNonInfectedLeaders) j];
    [size(+(minDistNonInfectedLeaders==1)) j];
    [size(+(minDistInfectedLeaders==1)) j];
    [size(+(minDistNonInfectedLeaders==1)) j];
    
    netstats(j).neighborOfInfected = ((+(minDistInfectedLeaders==1) - +(minDistNonInfectedLeaders==1))>0);
    netstats(j).neighborOfNonInfected = ((+(minDistInfectedLeaders==1) - +(minDistNonInfectedLeaders==1))<0);
    
    netstats(j).degree = sum(X,2);
    netstats(j).num_edges = sum(X(:));
    netstats(j).leaderneighborhood = ((leaders*ones(1,N) + ones(N,1)*leaders')>0);
    netstats(j).num_leaders = sum(leaders);
    netstats(j).num_leader_edges = (leaders'*X*leaders);
else % If we've already computed this once, let's not repeat ourselves.
    minDistFromLeaders = netstats(j).minDistFromLeaders;
    avgDistFromLeaders = netstats(j).avgDistFromLeaders;
    minDistInfectedLeaders = netstats(j).minDistInfectedLeaders;
    minDistNonInfectedLeaders = netstats(j).minDistNonInfectedLeaders;
end


switch version
    case 1
        % 1. Fraction of nodes that have no taking neighbors but are takers
        % themselves.
        infectedNeighbors = sum((ones(N,1)*infected').*X, 2); % number of infected neighbors
        
        if sum(infectedNeighbors==0 & netstats(j).degree>0)>0
            stats(1) = sum((infectedNeighbors==0 & infected==1 & netstats(j).degree>0))/sum(infectedNeighbors==0 & netstats(j).degree>0);
        elseif sum(infectedNeighbors==0 & netstats(j).degree>0)==0
            stats(1) = 0;
        end
        
        % 2. Fraction of individuals that are infected in the neighborhood of infected leaders stats(1)=0;
        if sum(netstats(j).neighborOfInfected)>0
            stats(2) = sum(infected.*netstats(j).neighborOfInfected)/sum(netstats(j).neighborOfInfected);
        else
            stats(2)=0;
        end
        
        % 3. Fraction of individuals that are infected in the neighborhood of non-infected leaders
        if sum(netstats(j).neighborOfNonInfected)>0
            stats(3) = sum(infected.*netstats(j).neighborOfNonInfected)/sum(netstats(j).neighborOfNonInfected);
        else
            stats(3)=0;
        end
        
        % 4. Covariance of individuals taking with share of neighbors taking
        NonHermits = (netstats(j).degree>0);
        ShareofTakingNeighbors = infectedNeighbors(logical(NonHermits))./netstats(j).degree(logical(NonHermits));
        NonHermitTakers = infected(logical(NonHermits));
        stats(4) = sum(NonHermitTakers.*ShareofTakingNeighbors)/sum(NonHermits);
        
        % 5. Covariance of individuals taking with share of second neighbors taking
        % taking
        infectedSecond = sum(Sec.*(infected*ones(1,N)),2);
        ShareofSecond = infectedSecond(logical(NonHermits))./netstats(j).degree(logical(NonHermits));
        stats(5) = sum(NonHermitTakers.*ShareofSecond)/sum(NonHermits);
        
        
    case 2
        % 1. Fraction of nodes that have no taking neighbors but are takers
        % themselves.
        infectedNeighbors = sum((ones(N,1)*infected').*X, 2); % number of infected neighbors
        
        if sum(infectedNeighbors==0 & netstats(j).degree>0)>0
            stats(1) = sum((infectedNeighbors==0 & infected==1 & netstats(j).degree>0))/sum(infectedNeighbors==0 & netstats(j).degree>0);
        elseif sum(infectedNeighbors==0 & netstats(j).degree>0)==0
            stats(1) = 0;
        end
        
        % 2. Covariance of individuals taking with share of neighbors taking
        NonHermits = (netstats(j).degree>0);
        ShareofTakingNeighbors = infectedNeighbors(logical(NonHermits))./netstats(j).degree(logical(NonHermits));
        NonHermitTakers = infected(logical(NonHermits));
        stats(2) = sum(NonHermitTakers.*ShareofTakingNeighbors)/sum(NonHermits);
        
        % 3. Covariance of individuals taking with share of second neighbors taking
        % taking
        infectedSecond = sum(Sec.*(infected*ones(1,N)),2);
        ShareofSecond = infectedSecond(logical(NonHermits))./netstats(j).degree(logical(NonHermits));
        stats(3) = sum(NonHermitTakers.*ShareofSecond)/sum(NonHermits);
        
        
    case 3
        % same as case 2, but purged ofleader injection points.
        leaderTrue = ((leaders > 0)); % a variable that denotes whether a node is either a leader
        
        % 1. Fraction of nodes that have no taking neighbors but are takers themselves.
        infectedNeighbors = sum((ones(N,1)*infected').*X, 2); % number of infected neighbors
        
        if sum(infectedNeighbors==0 & netstats(j).degree>0)>0
            stats(1) = sum((infectedNeighbors==0 & leaderTrue==0 & infected==1 & netstats(j).degree>0))/sum(infectedNeighbors==0 & leaderTrue==0 & netstats(j).degree>0);
        elseif sum(infectedNeighbors==0 & leaderTrue==0 & netstats(j).degree>0)==0
            stats(1) = 0;
        end
        
        
        % 2. Covariance of individuals taking with share of neighbors taking
        NonHermits = (netstats(j).degree>0);
        NonHermitsNonLeaders = ( NonHermits & (1-leaderTrue) ); % not isolates, not leaders
        
        ShareofTakingNeighbors = infectedNeighbors(logical(NonHermitsNonLeaders))./netstats(j).degree(logical(NonHermitsNonLeaders));
        NonHermitTakers = infected(logical(NonHermitsNonLeaders));
        stats(2) = sum(NonHermitTakers.*ShareofTakingNeighbors)/sum(NonHermitsNonLeaders);
        
        
        % 3. Covariance of individuals taking with share of second neighbors taking taking
        infectedSecond = sum(Sec.*(infected*ones(1,N)),2);
        ShareofSecond = infectedSecond(logical(NonHermitsNonLeaders))./netstats(j).degree(logical(NonHermitsNonLeaders));
        stats(3) = sum(NonHermitTakers.*ShareofSecond)/sum(NonHermitsNonLeaders);
        
    case 4
        % same as case 3, but purged of ALL leader nodes.
        leaderTrue = ((leaders > 0)); % a variable that denotes whether a node is either a leader
        
        % 1. Fraction of nodes that have no taking neighbors but are takers themselves.
        infectedNeighbors = sum((ones(N,1)*infected').*X*(1-leaderTrue), 2); % number of infected neighbors
        
        if sum(infectedNeighbors==0 & netstats(j).degree>0)>0
            stats(1) = sum((infectedNeighbors==0 & leaderTrue==0 & infected==1 & netstats(j).degree>0))/sum(infectedNeighbors==0 & leaderTrue==0 & netstats(j).degree>0);
        elseif sum(infectedNeighbors==0 & leaderTrue==0 & netstats(j).degree>0)==0
            stats(1) = 0;
        end
        
        
        % 2. Covariance of individuals taking with share of neighbors taking
        NonHermits = (netstats(j).degree>0);
        NonHermitsNonLeaders = ( NonHermits & (1-leaderTrue) ); % not isolates, not leaders
        
        ShareofTakingNeighbors = infectedNeighbors(logical(NonHermitsNonLeaders))./netstats(j).degree(logical(NonHermitsNonLeaders));
        NonHermitTakers = infected(logical(NonHermitsNonLeaders));
        stats(2) = sum(NonHermitTakers.*ShareofTakingNeighbors)/sum(NonHermitsNonLeaders);
        
        
        % 3. Covariance of individuals taking with share of second neighbors taking taking
        infectedSecond = sum(Sec.*(infected*ones(1,N))*(1-leaderTrue),2);
        ShareofSecond = infectedSecond(logical(NonHermitsNonLeaders))./netstats(j).degree(logical(NonHermitsNonLeaders));
        stats(3) = sum(NonHermitTakers.*ShareofSecond)/sum(NonHermitsNonLeaders);
        
        
end


return;
