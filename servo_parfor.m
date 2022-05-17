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
%
% SRH Notes on Weds Jan 13 2021: servo_parfor.m is SRH modification of
% RKP's servo.m code to do final algorithm runs in parallel. This file
% creates .mat output files that can then be used as inputs to the
% PlotQuantilesOracleName.m plotting code.
%==========================================================================

%==========================================================================
function servo()
    
    for numdim = 4:4 %1:2 %4
    for numalgs = 1:2
    
    %%%%%%%%%%%%%%%%%%%problem and algorithm parameters%%%%%%%%%%%%%%%%%
    %OracleName = 'BusScheduling';
    %OracleName = 'DiscreteQuadratic';
    OracleName = 'DynamicNews';
    if numalgs == 1
        AlgorithmName = 'ADALINE';
    elseif numalgs == 2
        AlgorithmName = 'RSPLINE';
    end
    
    %%
    if strcmp(OracleName,'BusScheduling') == 1
        if numdim == 1
            d = 9;                                 % numberof dimensions
            maxbudget = 20000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 2
            d = 20;
            maxbudget = 30000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 3
            d = 50;
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 4
            d = 100;
            maxbudget = 50000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        end
    end
    
    %%
    if strcmp(OracleName,'DiscreteQuadratic') == 1
        if numdim == 1
            d = 25;                                 % numberof dimensions
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 2
            d = 50;
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 3
            d = 100;
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 4
            d = 200;
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        end
    end
    
    %%
    if strcmp(OracleName,'DynamicNews') == 1
        if numdim == 1
            d = 5;                                  % number of dimensions
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 2
            d = 15;
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 3
            d = 25;
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        elseif numdim == 4
            d = 50;
            maxbudget = 40000;                      % maximum budget for each run
            fprintf('=============== STARTING %s d = %d ===============\n',AlgorithmName,d);
        end
    end
    
    %%%%%%%%%%%%%%plot and run parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(OracleName,'BusScheduling') == 1
        x_0 = ones(d,1);
    end
    if strcmp(OracleName,'DiscreteQuadratic') == 1
        x_0= 10 * ones(d,1);              % initial guess
    end
    if strcmp(OracleName,'DynamicNews') == 1
        x_0=  50 * ones(d,1);             % initial guess
    end
    reps = 1000;                            % number of macro reps of the algorithm
    tickpts =0:.02:1;                      % data collection pts
    %iseed_run = [360974; 84995; 89348; 578693; 35646; 7488595];                   % initial seeds for servo                              
    plotxs = round(maxbudget * tickpts);
    plotpaths = zeros(size(tickpts,2),reps+1);
    plotpaths(:,1) = plotxs';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% initialize
    %nextseed = iseed_run;

    %% run the algorithm a bunch of times and collect output
    wholerep = tic;
    parfor (j = 1:reps,12) %(j = 1:reps,72) %second number is max number of cores
        %%FOR USE WITH PARFOR COMMAND: GENERATE SEEDS RANDOMLY 
            u = rand(1,3);
            nextseed = round(10000000 * u' );
        %% obtain the initial pt
        % execute the algorithm and get the algorithmic path
        % alg_path is a m x 2 matrix; the first column lists the comp work, next d columns 
        % the returned solution, and the last column the optimality gap
        innerrep = tic;
        if ( AlgorithmName == 'ADALINE' )
            %fprintf('%s %d %s %s\n',"ADALINE starting replication ",j,"on ", OracleName);
            [alg_path] = adaline(nextseed, x_0, maxbudget,OracleName);
        elseif (AlgorithmName == 'RSPLINE' )
            %fprintf('%s %d %s %s\n',"RSPLINE starting replication",j,"on ", OracleName);
            [alg_path] = RSPLINE(nextseed, x_0, maxbudget,OracleName);
        end
        elapsed_time = toc(innerrep);
        fprintf('%s replication %3.0f took %2.0f:%2.0f \n',AlgorithmName,j, floor(elapsed_time/60), ceil(mod(elapsed_time,60)))
        %% find the discrete quantile (without interpolation;
        % determine the solutions that would have been returned by the
        % algorithm if the algorithm was stopped at plotxs
        [stoppedptsindex,~] = findinv(alg_path(:,1),plotxs');
        %% nextplotpath has 2 + d + 1 columns; the first two correspond to
        % request work and actually done work at stopping, the next d
        % crrespond to returned solution and the last the optimality gap
        nextplotpath = [plotxs',alg_path(stoppedptsindex,[1,(d+2)])];
        %% plotpaths a particular replication's data listed in a column; the
        % first column correspondsto the request work stop
        plotpaths(:,j+1) = nextplotpath(:,3);
    end
    elapsed_rep = toc(wholerep);
    fprintf('Whole Rep Time %3.1f\n',elapsed_rep/60);
    
    %%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% compute the quantile curves
    qtpts = [0.1,0.15,0.25,0.5,0.75,0.85,0.9];
    % quantilecurve is a matrix with the first column corresponding to the
    % work points and each of the rest of the columns corresponding to a
    % requested quantile
    quantilecurve = zeros(size(tickpts,2),size(qtpts,2)+1);
    quantilecurve(:,1) = plotxs';
    for w_index = 1 : size(tickpts,2)
        quantiles_w = quantile(plotpaths(w_index,2:reps+1),qtpts);
        quantilecurve(w_index,2:size(qtpts,2)+1) = quantiles_w;
    end
    
    dim = num2str(d);
    strreps = num2str(reps);
    myfilename=strcat('Data',OracleName,AlgorithmName,dim,'reps',strreps,'.mat');
    save(myfilename,'quantilecurve','qtpts');

    
    %% plot the quantile curves
%     plot(quantilecurve(:,1),quantilecurve(:,2));
%     xlim([0,maxbudget]);
%     ylim([-1,20]);
%     xlabel('Number of Oracle Calls');
%     ylabel('Quantile of True Optimality Gap %')
%     hold;
%     for j = 2 : size(qtpts,2)
%         plot(quantilecurve(:,1),quantilecurve(:,j+1));
%     end
    
    end
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

