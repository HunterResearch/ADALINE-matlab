function [seed,xbar,x2bar,feasibilityFLAG] = simulate(iseed,ssize,loc,feascheckonlyFLAG,OracleName)
    global pts_visited
    % if simulator is being called only to check feasibility, then do it
    % and quit
    seed = iseed; 
    xbar = nan;
    x2bar = nan;
    
    if ( feascheckonlyFLAG == 1 )
        [~,~,feasibilityFLAG] = call_oracle(iseed,loc,feascheckonlyFLAG,OracleName);
        return;
    end
    % find the initial seed to use; this is non trivial because we want to
    % use common random numbers
    % first find the location index of the point
    d = size(loc,1);
    if (isempty(pts_visited) == 0)
        [~,locindex,~] = intersect(pts_visited(:,1:d),loc','rows');
    else
        locindex = [];
    end
    if ( isempty(locindex) == 0 )
        % if the point exists, find the no of samples already observed 
        alreadysampled = pts_visited(locindex,(d+1));
        seed = iseed;
        % find the seed by running a uniform generator "alreadysampled"
        % times
        for k = 1 : alreadysampled
            [seed,~] = mrg32k3a(seed);
            [seed,~] = mrg32k3a(seed);
        end
    else
        seed = iseed;
    end
    % initialize
    sum = 0;
    sum2 = 0;
    % execute the oracle 
    for j = 1 : ssize
        feascheckonlyFLAG =0;
        [seed,nexty,feasibilityFLAG] = call_oracle(seed,loc,feascheckonlyFLAG,OracleName);
        sum = sum + nexty;
        sum2 = sum2 + nexty^2;
    end
    xbar = sum / ssize;
    x2bar = sum2 / ssize;
end