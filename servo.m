%==========================================================================
%                       ADALINE
%==========================================================================
% LAST UPDATE
%               Feb. 16 2020
% PROGRAMMER
%               Raghu Pasupathy
% REFERENCE		
%               Hunter, Pasupathy, Raghavan, and Taaffe (2019)
%		  		An Algorithm for Local Stochastic Optimization Over 
%               Integer Variables.
%==========================================================================
%
%
%==========================================================================
function servo()
    global OracleName
    
    %%%%%%%%%%%%%%%%%%%problem and algorithm parameters%%%%%%%%%%%%%%%%%
    %OracleName = 'BusScheduling';
    %OracleName = 'DiscreteQuadratic';
    OracleName = 'DynamicNews';
    AlgorithmName = 'adaline';
    %AlgorithmName = 'RSPLINE';
    d = 5;                                 % numberof dimensions 
    
    %%%%%%%%%%%%%%plot and run parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x_0=  50 * ones(d,1);                    % initial guess
    reps = 12;                              % number of macro reps of the algorithm
    maxbudget = 40000;                      % maximum budget for each run
    tickpts = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];    % data collection pts
    iseed_run = [8490106; 2973185; 3620756];                   % initial seeds for the solver                             
    plotxs = round(maxbudget * tickpts);
    plotpaths = zeros(size(tickpts,2),reps+1);
    plotpaths(:,1) = plotxs';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initialize
    %nextseed = iseed_run;

    % run the algorithm a bunch of times and collect output
    parfor j = 1 : reps
        u = rand(1,3);
        %u = [0.9303, 0.7517, 0.2979];
        nextseed = round(10000000 * u' );
        % obtain the initial pt
        % execute the algorithm and get the algorithmic path
        % alg_path is a m x 2 matrix; the first column lists the comp work, next d columns 
        % the returned solution, and the last column the optimality gap
        if ( AlgorithmName == 'adaline' )
            fprintf('%s %d %s %s\n',"Adaline starting replication ",j,"on ", OracleName);
            [alg_path] = adaline(nextseed, x_0, maxbudget,OracleName);
        elseif (AlgorithmName == 'RSPLINE' )
            fprintf('%s %d %s %s\n',"RSPLINE starting replication",j,"on ", OracleName);
            [alg_path] = RSPLINE(nextseed, x_0, maxbudget,OracleName);
        end
        
        % find the discrete quantile (without interpolation;
        % determine the solutions that would have been returned by the
        % algorithm if the algorithm was stopped at plotxs
        [stoppedptsindex,~] = findinv(alg_path(:,1),plotxs');
        % nextplotpath has 2 + d + 1 columns; the first two correspond to
        % request work and actually done work at stopping, the next d
        % crrespond to returned solution and the last the optimality gap
        nextplotpath = [plotxs',alg_path(stoppedptsindex,[1,(d+2)])];
        % plotpaths a particular replication's data listed in a column; the
        % first column correspondsto the request work stop
        plotpaths(:,j+1) = nextplotpath(:,3);
    end
    %%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the quantile curves
    qtpts = [0.25,0.5,0.75, 0.9];
    % quantilecurve is a matrix with the first column corresponding to the
    % work points and each of the rest of the columns corresponding to a
    % requested quantile
    quantilecurve = zeros(size(tickpts,2),size(qtpts,2)+1);
    quantilecurve(:,1) = plotxs';
    for w_index = 1 : size(tickpts,2)
        quantiles_w = quantile(plotpaths(w_index,2:reps+1),qtpts);
        quantilecurve(w_index,2:size(qtpts,2)+1) = quantiles_w;
    end
    
    % plot the quantile curves
    plot(quantilecurve(:,1),quantilecurve(:,2));
    xlim([0,maxbudget]);
    ylim([0,3]);
    xlabel('Number of Oracle Calls');
    ylabel('Quantile of True Optimality Gap %')
    hold;
    for j = 2 : size(qtpts,2)
        plot(quantilecurve(:,1),quantilecurve(:,j+1));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [locfinv,finv] = findinv(x,target)
% this function takes each element y intarget and finds the largest value in
% x that is smaller than or equal to y, that is, finds Finv(y)
    no_ys = size(target,1);
    no_xs = size(x,1);
    finv = -inf * ones(no_ys,1);
    locfinv = zeros(no_ys,1);
    target = sort(target);
    x = sort(x);
    j_k = 1;
    for j =  1: no_ys
        for k =  j_k : no_xs
            if ( x(k) <= target(j) && x(k) > finv(j) )
                finv(j) = x(k);
                locfinv(j) = k;
                j_knext = k;
            end
        end
        j_k = j_knext;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iseed,u16807d]=u16807d( iseed )
%..........................................................................
%     bruce schmeiser     january 1989.                                   .
%     a linear congruential pseudorandom number generator                 .
%       using constant 16807 and modulus (2**31)-1.                       .
%              iseed = iseed*16807 (mod 2^31 -1)                          .
%     for implementations that don't require double precision, see        .
%       s.k. park and k.w. miller, "random numbers generators: good       .
%         ones are hard to find," cacm, 31, 10 (October 1988), 1192-1201. .
%     in correct implementations, starting with a seed value of 1 will    .
%       result in a seed value of 1043618065 on call number 10001.        .
%..........................................................................
%     input:  iseed.   integer.                                           . 
%                        chosen from [1,2147483646] on the first call.    .
%                        thereafter, the value returned from the last call.
%     output: iseed.   integer.                                           .
%                        to be used in the next call.                     .
%     output: u16807d.  real.                                             .
%                        a pseudorandom number in (0,1).                  .
%..........................................................................
    %iseed=0;     
    u16807d=0;
    while (u16807d<=0 || u16807d>=1)
        iseed = mod (iseed * 16807,2147483647);
        u16807d = iseed / 2147483648;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
            
