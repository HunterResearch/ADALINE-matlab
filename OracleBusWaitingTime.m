function [nextseed, y, feasibilityFLAG]= OracleBusWaitingTime(nextseed, x, feascheckonlyFLAG)
    % this is the oracle for the bus scheduling problem
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  input:
    %      param = 
    %         param(1) = problem ID, 2 (not used!)
    %         param(2) = problem dimension
    %         param(3) = nseeds = 1
    %         param(4) = nsecMeas = 0
    %         param(5) = lambda, arrival rate of passengers
    %         param(6) = gamma, final departure time of bus
    %      x =
    %         matrix of size 1 X id 
    %         integer solution [x(1), x(2), ..., x(id)]
    %      m = sample size
    %      iseed = 
    %         matrix of size 1 X nseeds
    %         [iseed(1), ..., iseed(nseeds)]
    %  output:  
    %	   flag1 =      0 implies that the parameters are feasible, otherwise 
    %                   infeasible
    %      flag2 =      0 implies that x is feasible, otherwise infeasible
    %      fn    =      ybar (defined only if flag1 = 0 and flag2 = 0)	
    %      constraint = matrix (size = 1 X nsecMeas) of estimates of 
    %                   constraint functions 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initialize and set problem parameters 
        m=1;
        iseed = nextseed;
        id     = size(x,1);
        lambda = 10;
        gamma  = 100;
        % this oracle requires an input that is a row vector!
        x = x';
        y = nan;

        %% model parameters feasibility check
        flag1 = 0;
        if lambda < 0 || gamma < 0 || id~=length(x) 
            flag1 = 1;
        end

        % check the feasibility of x
        feasibilityFLAG = 1;
        if ( min(x) < 0 || max(x) > gamma ) 
            feasibilityFLAG = 0;
            y = inf;
        end
        % if you are just checking feasibility, get out
        if ( feascheckonlyFLAG == 1 )
            return;
        end

        %decision-variable feasibility check
        flag2 = 0;
        if (sum(x < 0) > 0 || sum(x>gamma) > 0)
            flag2=1;
        end

        if ( flag1 == 1 || flag2 == 1)
            % infeasible and so return infinity
            y = inf;
            return
        end

        sum1 = 0;
        sum2 = 0;
        for i=1:m
            timesum = 0;
            tarrive = 0;
            timebus = 0;

            %total wait of all passengers on day i
            while (1)
                [iseed, u]= mrg32k3a(iseed);
                if(1- u > 0)
                    tarrive = tarrive + ((-log(1.0 - u))/lambda);
                end
                if(tarrive > gamma) 
                    break
                end
                if(tarrive > timebus)  %compute time of the next bus
                    timebus = gamma;
                    for j=1:id	
                        if(x(j) >= tarrive && x(j) < timebus)
                            timebus=x(j);
                        end
                    end
                end
                timesum = timesum + (timebus - tarrive);
            end
            sum1=sum1+timesum;
            sum2=sum2+timesum^2;
        end
        y=sum1/m;
        nextseed = iseed;
end