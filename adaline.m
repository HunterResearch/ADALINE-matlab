%==========================================================================
%                       ADALINE
%==========================================================================
% LAST UPDATE
%               Feb. 15 2020
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
function [alg_path] = adaline(iseed,x_0,budget, OracleName)
    % read inputs.
    d = size(x_0,1);
    b = min(sqrt(d),5);                  % max number of line searches during each iteration
    
    % use iseed(1) to create y_seed
    seed = iseed(1,1); y_seed = zeros(6,1);
    for j = 1 : 6
        [seed, u] = u16807d( seed );        
        y_seed(j,1) = floor(10^6 * u );
    end
    
    % use iseed(2) to create daseed
    seed = iseed(2,1); daseed = zeros(6,1);
    for j = 1 : 6
        [seed, u] = u16807d( seed );        
        daseed(j,1) = floor(10^6 * u );
    end
    
    % use iseed(3) to create ptseed
    seed = iseed(3,1); pt_seed = zeros(6,1);
    for j = 1 : 6
        [seed, u] = u16807d( seed );        
        pt_seed(j,1) = floor(10^6 * u );
    end
    
    alg_path = [[0,x_0'],find_optgap(x_0,OracleName)];
    
    % initialize
    iteration=1;
    X_kminusone = x_0;
    Xtilde_k = X_kminusone;
    workdone = 0;
    neMtilde_k = 2;
    necalls = 0;
    liMtilde_k = 0;
   

    % pts_visited is a global database of visited points; it is legacy and
    % it plays no role here
    global pts_visited;
    global budgetexceed_FLAG
    pts_visited= [];
    
    % continue until simulation budget is exhausted
    while ( workdone < budget )
        % update the escort sequence
        lambda_k = ceil(max(2 * log(iteration),3));
        budgetrem = budget - workdone;

        % advance the start seed if necessary to perform independent sampling across NEs 
        advanceby = max(neMtilde_k,liMtilde_k);
        j = 0;
        while ( j < advanceby ) 
            [y_seed, ~] = mrg32k3a(y_seed);
            [daseed, ~] = mrg32k3a(daseed);
            j = j + 1;
        end
        
        fromneFLAG = 0;
        
        if iteration == 1
            [~,nextxbar,~,~] = simulate(y_seed,neMtilde_k,X_kminusone,0,OracleName);
            FXtilde_k = nextxbar;
            Xtilde_k = X_kminusone;
        end
 
        % start with NE from the second iteration onward
        if iteration > 1
            [pt_seed,~, FXtilde_k, Xtilde_k,neMtilde_k,calls] = NE(pt_seed, y_seed, X_kminusone, iteration, lambda_k, budgetrem,OracleName);
            fromneFLAG = 1;
            necalls = calls;
            workdone = workdone + calls;
        
            % record progress
            f_x = find_optgap(Xtilde_k,OracleName);
            alg_path = [alg_path; [workdone, Xtilde_k',f_x]];
            
        end
            
        M_k = max(neMtilde_k,lambda_k);
        
        % this DA call is substantive only in the first iteration since NE is
        % not performed in the first iteration. (DA is called within LI.)
        [~,~,dhat_k, Xtilde_k, daMtilde_k,dacalls,timetoperformNEFLAG] = DA(daseed,y_seed,X_kminusone, Xtilde_k,FXtilde_k,M_k,lambda_k,iteration, budget - workdone, fromneFLAG,OracleName);
        workdone = workdone + dacalls;
        
        % record progress
        f_x = find_optgap(Xtilde_k,OracleName);
        alg_path = [alg_path; [workdone, Xtilde_k',f_x]];
        
        % perform LI
        [~,~,X_kminusone,calls,~,totdacalls,~, N_k,liMtilde_k] = LI(daseed,y_seed, FXtilde_k, Xtilde_k, dhat_k, daMtilde_k, lambda_k, b, budget - workdone, iteration,OracleName,timetoperformNEFLAG);
        workdone = workdone + calls;
        
        % record progress
        f_x = find_optgap(X_kminusone,OracleName);
        alg_path = [alg_path; [workdone, X_kminusone',f_x]];

        iteration = iteration + 1;
        
        % REPORT PROGRESS
        %l_inf = max(abs(X_kminusone));
        %fprintf('%4d %10.3f %8d %8d %4d %8d %4d %8d \n',iteration-1, f_x, max(X_kminusone) - min(X_kminusone), necalls, neMtilde_k, N_k, totdacalls, workdone);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pt_seed, y_seed, FXtilde_k, Xtilde_k,Mtilde_k,calls] = NE(pt_seed, y_seed, X_kminusone, iteration, lambda_k, nebudget,OracleName)
    global pts_visited
    global budgetexceed_FLAG
    
    % algorithmc constant; decides when to quit NE
    t_critical = 1.65;
    
    %initialize
    calls = 0;
    d = size(X_kminusone,1);
    ssizemin = 1;
    FXtilde_k = inf;
    Mtilde_k = 0;
    Xtilde_k = X_kminusone;
    budgetexceed_FLAG = 0;
    pts_visited = [];
    m = 0;
    
    % create a data structure for storing neighbor data
    % the jth object in nbrdata is a n x 2 column of data; the first column corresponds to data observed in the 
    % j+ direction and the second to the j- direction
    for j = 1 : d
        nbrdata{j} = [];
    end
    
    better_nbr_FLAG = 0;            % better neighbor flag is false
    
    % delnbdN1 is the deleted N1 neighborhood of X_kminusone; delnbdN1 is a
    % matrix with at most 2d columns and exactly d rows; 
    [delnbdN1] = find_delnbdN1(X_kminusone,OracleName);     
    no_delnbdN1 = size(delnbdN1,2);
    
    A = []; Ainfo=[]; 
    center_data = [];                   % all data at the center
    pseudo_der = zeros( d, 4 );         % pseudo derivatives plus and minus in columns 1,3; standard errors in 2,4
    
    
    % create a pmf supported on the neighborhood; mu_klist is a matrix
    % having no_delnbdN1 rows and (d+1) columns; the first d columns
    % correspond to the coordinates and the last column to the probability
    % mass
    mu_klist(:,1:d) = delnbdN1';
    mu_klist(:,(d+1)) = 1 / no_delnbdN1;

    % now try to find a better nbr
    while ( ( ( better_nbr_FLAG == 0 ) | ( Mtilde_k < lambda_k ) ) & ( calls < nebudget ) )
        % sample a point W from the support of A
        [pt_seed,W] = gendelnbd(pt_seed,mu_klist);
        % sample an observation from the pair corresponding to W
        feasonlycheckFLAG = 0;
        deviation = round(abs(W - X_kminusone));
        [~,dirindex] = max(deviation); 
        % simulate the pair
        [~,nextxplus,~,plusfeasFLAG] = simulate(y_seed,1,X_kminusone + deviation,feasonlycheckFLAG,OracleName);
        [~,nextxminus,~,minusfeasFLAG] = simulate(y_seed,1,X_kminusone - deviation,feasonlycheckFLAG,OracleName); 
        nbrdata{dirindex} = [nbrdata{dirindex};[nextxplus, nextxminus]];
        m = size(nbrdata{dirindex},1);
        
        % update the function values, global and local lists
        if (plusfeasFLAG == 1)
            [ssize, nextxbarplus,nextx2barplus] = updategloballist(X_kminusone + deviation,1,nextxplus,nextxplus^2,iteration);
            [A, Ainfo] = updatelocallist(A,Ainfo,X_kminusone + deviation,ssizemin,ssize,nextxbarplus,nextx2barplus,iteration);
        end
        
        if (minusfeasFLAG == 1)
            [ssize, nextxbarminus, nextx2barminus] = updategloballist(X_kminusone - deviation,1,nextxminus,nextxminus^2,iteration);
            [A, Ainfo] = updatelocallist(A,Ainfo,X_kminusone - deviation,ssizemin,ssize,nextxbarminus,nextx2barminus,iteration);
        end
        
        % simulate the center if needed
        Deltastar = m - size(center_data,1);
        if ( Deltastar > 1 )
            frpintf("Something seems to be wrong here ... CRN not in effect!");
        end
        if ( Deltastar > 0 )
            [~,nextx,~,~] = simulate(y_seed,1,X_kminusone,feasonlycheckFLAG,OracleName);
            center_data = [center_data;nextx];
            [~, ~, ~] = updategloballist(X_kminusone,1,nextx,nextx^2,iteration);
        end
        
        calls = calls + ( plusfeasFLAG + minusfeasFLAG + max(Deltastar,0) );
        
        % update the pseudo derivative estimates and their standard errors; use partial CRN 
        center_ssize = size( center_data, 1);
        F_Xkminusone = sum( center_data( : , 1 ) ) / center_ssize;
        pseudo_der_plus = (1/2) * ( ( sum(nbrdata{dirindex}(1:m,1))  / m ) - ( sum(nbrdata{dirindex}(1:m,2))  / m ) );
        pseudo_der_minus = -pseudo_der_plus;
        devplus = ( 1/2 ) * ( nbrdata{dirindex}(m,1) - nbrdata{dirindex}(m,2) );
        devminus = -devplus;
        pseudo_der(dirindex,1) =  pseudo_der_plus;
        pseudo_der(dirindex,2) = ( pseudo_der(dirindex,2) * ( m - 1 )  + devplus^2 ) / m;
        pseudo_der(dirindex,3) =  pseudo_der_minus;
        pseudo_der(dirindex,4) = ( pseudo_der(dirindex,4) * ( m - 1 )  + devminus^2 ) / m;      
        
        % identify the candidate best point
        if (isempty(Ainfo) == 0 )
            [FXtilde_k,argminestindex] = min(Ainfo(:,d+2));
            Xtilde_k = (Ainfo(argminestindex,(1:d)))';
            Mtilde_k = Ainfo(argminestindex,( d + 1 ));
        end
        
        % update the directional derivatives along the 2d directions
        dirindex_ssize = zeros(d,1);
        varvectorplus = inf(d,1);
        varvectorminus = inf(d,1);
        for dirindex = 1 : d
            dirindex_ssize(dirindex) = size( nbrdata{dirindex},1);
            % calculate the variance
            if dirindex_ssize(dirindex) > 1
                varvectorplus(dirindex) = ( pseudo_der(dirindex,2) - ( pseudo_der(dirindex,1) .* pseudo_der(dirindex,1) ) ) ./ dirindex_ssize(dirindex) ;
                varvectorminus(dirindex) = ( pseudo_der(dirindex,4) - ( pseudo_der(dirindex,3) .* pseudo_der(dirindex,3) ) ) ./ dirindex_ssize(dirindex) ;
            end
        end
        Tscoresplus = max(-pseudo_der(:,1), 0) ./ ( sqrt( varvectorplus ) );
        Tscoresminus = max(-pseudo_der(:,3), 0) ./ ( sqrt( varvectorminus ) );

        % clean up the nans
        for dirindex = 1 : d
            %if varvectorplus(dirindex) == 0
            %    Tscoresplus(dirindex) = 0;
            %end
            if isnan(Tscoresplus(dirindex))
                Tscoresplus(dirindex) = 0;
            end
            %if varvectorminus(dirindex) == 0
            %    Tscoresminus(dirindex) = 0;
            %end
            if isnan(Tscoresminus(dirindex))
                Tscoresminus(dirindex) = 0;
            end
        end

        weighTsplus = zeros(d,1);
        weighTsminus = zeros(d,1);
        % assign weight if not NaN
        for dirder = 1 : d
            if Tscoresplus(dirder) > 0 
                weighTsplus(dirder) = Tscoresplus(dirder);
            end
            if Tscoresminus(dirder) > 0 
                weighTsminus(dirder) = Tscoresminus(dirder);
            end
        end
        
        %%%%%%%%%%%%% legacy code commented %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Xbardiffplus = eye(d); Xbardiffminus = -eye(d);
        %dhat_k = ( Xbardiffminus * weighTsminus ) + ( Xbardiffplus * weighTsplus );
        %if ( norm(dhat_k) == 0 )
        %    [maxplus,maxlocplus] = max(weighTsplus);
        %    [maxminus,maxlocminus] = max(weighTsminus);
        %    if maxplus > maxminus
        %        dhat_k = Xbardiffplus(:,maxlocplus);
        %    else
        %        dhat_k = Xbardiffminus(:,maxlocminus);
        %    end            
        %else
        %    dhat_k = dhat_k / norm( dhat_k );
        %end
        
        if ( max( max( weighTsplus ), max( weighTsminus ) ) > t_critical )
            better_nbr_FLAG = 1;
        else
            better_nbr_FLAG = 0;
        end
        
        % update pmf
        %X_kminusoneinfo = findinfo(X_kminusone);
        %if isnan(dhat_k)
        %    fprintf("hello");
        %end
        mu_klist = updatepmf(mu_klist,Ainfo,center_data);
    end    
    pts_visited= [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delnbdx] = find_delnbdN1(x,OracleName)
    % cycle through each of the coordinates and gather all the feasible
    % neighbors; delnbdx has d rows and as many columns as neighbors
    d = size(x,1);
    delnbdx = [];
    basisvectors = eye(d);
    feascheckonlyFLAG = 1;
    for coord = 1:d
        nextbasis =basisvectors(:,coord);
        [~,~,~,feasibilityFLAG] = simulate(1,1,x+nextbasis,feascheckonlyFLAG,OracleName);
        if (feasibilityFLAG == 1)
            delnbdx = [delnbdx,x+nextbasis];
        end
        [~,~,~,feasibilityFLAG] = simulate(1,1,x-nextbasis,feascheckonlyFLAG,OracleName);
        if (feasibilityFLAG == 1)
            delnbdx = [delnbdx,x-nextbasis];
        end  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [extract] = findinfo(M)
    global pts_visited
    extract = [];
    % M is a matrix with each column corresponding to a point; M has d rows
    % find the number of pts
    pts = size(M,2);
    if isempty(pts_visited) == 0
        d = size(pts_visited,2) - 4;
        % extract all the rows in the global list corresponding to M
        for j = 1 : pts
            x = M(:,j);
            [~,iglobal,~] = intersect(pts_visited(:,1:d),x','rows');
            if (isempty(iglobal) == 0)
                extract = [extract; pts_visited(iglobal,:)];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ssize_updated, xbar, x2bar] = updategloballist(x,ssize,xbar,x2bar,iteration)
    global pts_visited
    d = size(x,1);
    % find the row that should be updated/populated and tackle the first d cols
    % first check if the global list is empty
    if ( isempty(pts_visited) == 0)
        [~,iglobal,~] = intersect(pts_visited(:,1:d),x','rows');
        % global list is nonempty but there exists no point match
        if ( isempty(iglobal) == 1 )
            iglobal = size(pts_visited,1) + 1;
            pts_visited(iglobal, 1:d) = x;
            pts_visited(iglobal, (d+4)) = iteration;
        end        
    else
        iglobal = 1;
        pts_visited(iglobal,1:d) = x;
        pts_visited(iglobal,(d+4)) = iteration;
    end
    % now update the rest of the columns
    % obtain the existing values first
    temp = pts_visited(iglobal,(d+1):(d+3)); 
    ssize_old = temp(1); xbar_old = temp(2); x2bar_old = temp(3);
    ssize_updated = ssize_old + ssize;
    xbar = ( ( xbar_old * ssize_old ) + ( xbar * ssize ) ) / (ssize_old + ssize);
    x2bar = ( (x2bar_old * ssize_old) + ( x2bar * ssize ) ) / (ssize_old + ssize);
    pts_visited(iglobal,d+1) = ssize_updated;
    pts_visited(iglobal,d+2) = xbar;
    pts_visited(iglobal,d+3) = x2bar;
    pts_visited(iglobal,d+4) = iteration;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, Ainfo] = updatelocallist(A, Ainfo,x,ssizemin,ssize,xbar,x2bar,iteration)
    ilocal = [];
    d = size(x,1);
    % update A and Ainfo only if the sample size is larger than the minimum
    % for update
    if (ssize >= ssizemin )
        if ( isempty(Ainfo) == 0 )
            % find the row that should be updated
            [~,ilocal,~] = intersect(Ainfo(:,1:d),x','rows');
        end
        % if there exists a row already nothing needs be done to A; just update Ainfo
        % otherwise, add a row A and Ainfo 
        if ( isempty(ilocal) == 0 )
            % update the row with the index ilocal
            Ainfo(ilocal,d+1) = ssize;
            Ainfo(ilocal,d+2) = xbar;
            Ainfo(ilocal,d+3) = x2bar;
            Ainfo(ilocal,d+4) = iteration;
        elseif( ssize >= ssizemin )
            % add a pt (column) to A
            A = [A,x];
            % add a row to Ainfo
            Ainfo = [Ainfo; [x',ssize,xbar,x2bar,iteration]];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [seed,w] = gendelnbd(seed,pmf)
    % we want to generate a point from pmf
    % mu_klist is a matrix having (no_delnbd1+1) rows and (d+1) columns; 
    % the first d columns correspond to the coordinates and the last column 
    % to the probability mass
    % generate a point from pmf
    % We need to use the alias method here for better performance.
    % I am being lazy, not implementing the alias method
    d = size(pmf,2) - 1;
    [seed,u] = mrg32k3a( seed );
    cdf = cumsum(pmf(:,d+1));
    w_index = sum(cdf < u * cdf(end)) + 1;
    w = (pmf(w_index,1:d))';
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
function [pmf] = updatepmfold(pmf,Ainfo, F_k)
% this is a module that computes a subjective pmf representing the "best pt" on the deleted
% nbd of a given point; the pmf is computed based on a t-statistic analogue
% 
% it is assumed that Ainfo is an extract from the global visited matrix. Ainfo has
% exactly (d+4) columns; the first d columns correspond to the coordinates; 
% (d+1) corresponds to the sample size, (d+2) to function estimate, (d+3)
% to raw second mean, and (d+4) to the last iteration. Ainfo has as many rows
% as there are points in the delc=eted nbd with sample size at least ssizemin.
% bestpt, as the name suggests, is the vector with the
% coordinates of the best point in M
%
% pmf is a |del. nbd| x (d+1) matrix; each row corresponds to a point in the deleted
% nbd; the (d+1)th column corresponds to a probability.

    % algorithm constants
    p_0 = 0.5;      % the probability of an unobserved system being better
    
    d = size(Ainfo,2) - 4;
    no_delnbd = size(pmf,1);
    cardAinfo = size(Ainfo,1);
    
    % add five columns to temp; the columns d+1 through d+4 will contain
    % the information from Ainfo; the column d+5 will contain the
    % probability of sampling the point
    pmftemp = [pmf(:,1:d),zeros(no_delnbd,1), inf(no_delnbd,1), inf(no_delnbd,1), inf(no_delnbd,1), pmf(:,(d+1))];
    [~,observedindexespmftemp,observedindexesAinfo] = intersect(pmftemp(:,1:d),Ainfo(:,1:d),'rows');
    pmftemp(observedindexespmftemp,(d+1):(d+4)) = Ainfo(observedindexesAinfo,(d+1):(d+4));
    pmftemp(observedindexespmftemp, (d+5)) = ones(size(observedindexespmftemp,1),1);
    
    % find the minium estimate and the argmin index in pmftemp
    minest = min(pmftemp(:,d+2));
    argminestindexes = find(pmftemp(:,d+2)==minest); 
    
    % find the indexes of observed suboptimal and unobserved
    observedindexessuboptimalpmftemp = setdiff(observedindexespmftemp,argminestindexes);
    unobservedindexespmftemp = setdiff([1:no_delnbd]', observedindexespmftemp);
    
    % compute all the optimality gap estimates
    gaps = abs(pmftemp(observedindexessuboptimalpmftemp,d+2) - minest);
    
    % compute the scores and relprobs for the suboptimal points
    sigma2s = pmftemp(observedindexessuboptimalpmftemp,(d+3)) - ( pmftemp(observedindexessuboptimalpmftemp,(d+2)) ).^2;
    vars = sigma2s ./ pmftemp(observedindexessuboptimalpmftemp,(d+1));
    scores =  (1/2) * ( gaps.^2 ) ./ vars;
    pratios = scores .^(-1);
    relprobs = pratios / sum( pratios );
    
    % compute the variance at the best point
    sigma2best = pmftemp(argminestindexes,(d+3))   - ( pmftemp(argminestindexes,(d+2)) ).^2 ;
    varbest = sigma2best ./ pmftemp(argminestindexes,(d+1));
    varbesteq = sum( varbest .^ -1 ) .^ -1;
    hbesteq = sqrt( varbesteq * sum( ( relprobs.^2 ) ./ vars ) );
   
    
    % calculate the probability of sampling from the already observed group
    % (exploitation probability)
    E_N_betterobs = sum(max(sign(F_k - pmftemp(observedindexespmftemp,d+2)),0));
    E_N_betterunobs = p_0 * ( no_delnbd - cardAinfo );
    % prob_obs = ( cardAinfo / no_delnbd );
    if ( E_N_betterunobs == 0 )
        prob_obs = 1;
    else
        prob_obs = E_N_betterobs / ( E_N_betterobs + E_N_betterunobs );
        %prob_obs = 0;
    end
    
    
    
    % compute the probability of sampling from each in the best group given that you are
    % sampling from the already observed group
    groupprob_best_given_observed = hbesteq / ( 1 + hbesteq );
    indivprobs_best_given_observed = groupprob_best_given_observed * varbest / sum( varbest ); 
    
    % compute the probability of sampling from each in the suboptimal group given that you are
    % sampling from the already observed group
    indivprobs_suboptimal_given_observed = ( 1 - groupprob_best_given_observed ) * pratios / sum( pratios );
    
    % now re-assign the probabilities
    indivprobs_best = prob_obs * indivprobs_best_given_observed;
    indivprobs_suboptimal = prob_obs * indivprobs_suboptimal_given_observed;
    pmftemp(observedindexessuboptimalpmftemp,d+5) = indivprobs_suboptimal;
    pmftemp(argminestindexes,(d+5)) = indivprobs_best;
    pmftemp(unobservedindexespmftemp,d+5) = (1-prob_obs) / size(unobservedindexespmftemp,1);
    pmf = pmftemp(:,[1:d,(d+5)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pmf] = updatepmf(pmf,Ainfo, x0_data)
% this is a module that computes a subjective pmf representing the "best pt" on the deleted
% nbd of a given point; the pmf is computed based on a t-statistic analogue
% 
% it is assumed that Ainfo is an extract from the global visited matrix. Ainfo has
% exactly (d+4) columns; the first d columns correspond to the coordinates; 
% (d+1) corresponds to the sample size, (d+2) to function estimate, (d+3)
% to raw second mean, and (d+4) to the last iteration. Ainfo has as many rows
% as there are points in the delc=eted nbd with sample size at least ssizemin.
% bestpt, as the name suggests, is the vector with the
% coordinates of the best point in M
%
% pmf is a |del. nbd| x (d+1) matrix; each row corresponds to a point in the deleted
% nbd; the (d+1)th column corresponds to a probability.

    % algorithm constants
    p_0 = 0.5;      % the probability of an unobserved system being better
    
    d = size(Ainfo,2) - 4;
    no_delnbd = size(pmf,1);
    cardAinfo = size(Ainfo,1);
    
    % get the center information
    n = size(x0_data,1);
    F_k = sum( x0_data ) / n;
    sigmaF_k = max(sum( (x0_data .* x0_data) - ( sum(x0_data) / n )^2 ) / n,0);
    
    
    % add five columns to temp; the columns d+1 through d+4 will contain
    % the information from Ainfo; the column d+5 will contain the
    % probability of sampling the point
    pmftemp = [pmf(:,1:d),zeros(no_delnbd,1), inf(no_delnbd,1), inf(no_delnbd,1), inf(no_delnbd,1), pmf(:,(d+1))];
    [~,observedindexespmftemp,observedindexesAinfo] = intersect(pmftemp(:,1:d),Ainfo(:,1:d),'rows');
    pmftemp(observedindexespmftemp,(d+1):(d+4)) = Ainfo(observedindexesAinfo,(d+1):(d+4));
    pmftemp(observedindexespmftemp, (d+5)) = ones(size(observedindexespmftemp,1),1);
    
    % find the minium estimate and the argmin index in pmftemp
    %minest = min(pmftemp(:,d+2));
    %argminestindexes = find(pmftemp(:,d+2)==minest); 
    
    % find the indexes of observed suboptimal and unobserved
    %observedindexessuboptimalpmftemp = setdiff(observedindexespmftemp,argminestindexes);
    unobservedindexespmftemp = setdiff([1:no_delnbd]', observedindexespmftemp);
    
    % compute all the optimality gap estimates
    %gaps = abs(pmftemp(observedindexessuboptimalpmftemp,d+2) - minest);
    
    % compute the Tscores for all observed points
    % objdelta = ( sum(x0_data) ./ pmftemp(observedindexespmftemp,(d+1)) )  - pmftemp(observedindexespmftemp,d+2);
    objdelta = F_k  - pmftemp(observedindexespmftemp,d+2);
    % compute the sigma^2 estimate; negative values are set to zero
    sigma2s = max( pmftemp(observedindexespmftemp,(d+3)) - ( pmftemp(observedindexespmftemp,(d+2)) ).^2, 0 );
    vars = ( sigma2s + sigmaF_k ) ./ pmftemp(observedindexespmftemp,(d+1));
    relprobstemp = tcdf( objdelta ./ sqrt(vars), pmftemp(observedindexespmftemp,(d+1)) ); 
    relprobs = relprobstemp / sum(relprobstemp);
    
    % calculate the probability of sampling from the already observed group
    % (exploitation probability)
    E_N_betterobs = sum(max(sign(F_k - pmftemp(observedindexespmftemp,d+2)),0));
    E_N_betterunobs = p_0 * ( no_delnbd - cardAinfo );
    % prob_obs = ( cardAinfo / no_delnbd );
    if ( E_N_betterunobs == 0 )
        prob_obs = 1;
    else
        prob_obs = E_N_betterobs / ( E_N_betterobs + E_N_betterunobs );
        %prob_obs = 0;
    end
    
    % compute the probability of sampling from each in the best group given that you are
    % sampling from the already observed group
    pmftemp(unobservedindexespmftemp,d+5) = (1-prob_obs) / size(unobservedindexespmftemp,1);
    pmftemp(observedindexespmftemp,d+5) = prob_obs * relprobs;
    pmf = pmftemp(:,[1:d,(d+5)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [daseed, y_seed, dhat_k, better_nbr, Mtilde_k, calls,timetoperformNEFLAG] = DA(daseed, y_seed, X_k, Xtilde_k,F_k,Mtilde_k,lambda_k, iteration, dabudget, fromneFLAG,OracleName)
% this module returns a unit direction vector that is a descent direction in a probabilistic sense; it also 
% returns a sample size to perform a line search along that direction  
    global budgetexceed_FLAG;
    global pts_visited
    pts_visited = [];
    % initialize 
    calls = 0;
    dhat_k = Xtilde_k - X_k;
    timetoperformNEFLAG = 0;
    descentfoundFLAG = 0;
    d = size(X_k,1);
    better_nbr = Xtilde_k; 
    
    % generate a 2d x d matrix of orthogonal basis vectors stored along rows; 
    % the first d rows correspond to e_j and the
    % last d rows correspond to -e_{j - d} for j = d+1, d+2, ...
    basisrows = zeros( (2*d), d );
    seed = daseed;
    % first generate (d-1) rows
    for j = 1 : d
        % generate a random number 1,-1 
        flip = 1;
        [seed,u] = mrg32k3a(seed);
        if ( u <= 0.5 )
            flip = -1;
        end
        basisrows(j,j) = flip;
        basisrows(d+j,j) = -flip;
    end
    
    if fromneFLAG == 1
        return;
    end
    
    % construct the descent cone 
    j = 1; 
    nextpt = X_k;
    currentbest = X_k;
    currentbestfn = F_k;
    fullbody = [];
    feascheckonlyFLAG =0 ;
    [~,nextxbar,nextx2bar,~] = simulate(y_seed,Mtilde_k,nextpt,feascheckonlyFLAG,OracleName);
    calls = calls + Mtilde_k;
    oldxbar = nextxbar;
    oldx2bar = nextx2bar;
    fullbodycount = 1;
    dirders = []; dirders_se = [];
    descentcone = [];
    nearconstraintFLAG = 0;
    Tscores=-Inf;
    
    while ( ( fullbodycount < (d + 1) ) & j <= ( 2 * d ) & ( max(Tscores) <= inf | isnan(Tscores) ) )
        nextpt = nextpt + basisrows(j,1:d)';
        [~,nextxbartemp,nextx2bartemp,feasFLAG] = simulate(y_seed,Mtilde_k,nextpt,feascheckonlyFLAG,OracleName);
        if ( feasFLAG == 1 )
            nextxbar = nextxbartemp;
            nextx2bar = nextx2bartemp;
            fullbodycount = fullbodycount + 1;
            calls = calls + Mtilde_k;
            fullbody = [fullbody, basisrows(j,1:d)'];
            dirders = [dirders; (nextxbar - oldxbar)];
            descentcone = [descentcone; min(nextxbar - oldxbar,0)];
            dirders_se = [dirders_se; sqrt(( ( nextx2bar - nextxbar^2 )  + ( oldx2bar - oldxbar^2 ) ) / Mtilde_k)];
            Tscores= -descentcone ./ dirders_se;
            oldxbar = nextxbar;
            oldx2bar = nextx2bar;
            % update the best point if needed
            if nextxbar < currentbestfn
                currentbest = nextpt;
                currentbestfn = nextxbar;
            end
        else
            nearconstraintFLAG = 1;
        end
        better_nbr = currentbest;
        j = j + 1;
    end
    % update direction (if near a constraint go along descent cone,
    % otherwise negative gradient
    if min(dirders) == 0
        dhat_k = [];
    elseif fullbodycount == (d+1)
        dhat_k = -( fullbody * dirders ) / norm( fullbody * dirders );
    else
        %dhat_k = -( fullbody * ( descentcone ./ dirders_se ) ) / norm( fullbody * ( descentcone ./ dirders_se ) );
        dhat_k = -( fullbody * descentcone ) / norm( fullbody * descentcone );
    end
    %distbest = norm( Xtilde_k - currentbest );
    if isempty(dhat_k) == 1
        dhat_k = zeros(d,0);
        timetoperformNEFLAG = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [daseed,y_seed,X_k,calls,totlicalls,totdacalls,RAcalls, N_k, maxSampleSize] = LI(daseed, y_seed,FXtilde_k, Xtilde_k, dhat_k, Mtilde_k, lambda_k, b, LIbudget,iteration,OracleName, timetoperformNEFLAG)
    global budgetexceed_FLAG
    global pts_visited
    pts_visited = [];
    d = size(Xtilde_k,1);
    % initialize module constants 
    l0 = 3 * sqrt(d);
    clickbacktol = 1;
    % initialize
    X_k = Xtilde_k;
    M_k = Mtilde_k;
    F_k = FXtilde_k;
    Fnext = -Inf;
    step_factor = 2;
    maxSampleSize = 0;
    expandstepsMAX =inf;
    eta = 0.9;
        
    budgetexceed_FLAG = 0;
    N_k = 0;
    calls = 0;
    totlicalls = 0;
    totdacalls = 0;
    RAcalls = 0;
    
    if ( timetoperformNEFLAG == 1 )
        return;
    end
    
    % perform at most b line searches or until local min is found
    while ( (N_k < b) && timetoperformNEFLAG == 0 )
        % initialize step size
        l = l0;
        % perform the next line search starting at W_k and with samp size M_k
        % do the expansion phase
        expandFLAG = 1;
        expandsteps = 0;
        
        % simulate at the start point of line search
        [~,F_k,~,~] = simulate(y_seed,M_k,X_k,0,OracleName);
        
        while ( expandFLAG == 1 && expandsteps < expandstepsMAX)
            % identify a feasible Wnext
            projectfeasibleFLAG = 0; multiplier = l;
            while projectfeasibleFLAG == 0
                Wnext = round( X_k + ( multiplier * dhat_k ) );
                % check feasibility
                [~,~,~,projectfeasibleFLAG] = simulate(1,1,Wnext,1,OracleName);
                % warn
                %if ( projectfeasibleFLAG == 0 )
                %    fprintf("%s \n", "Warning: suggested point is infeasble ... backtracking");
                %end
                multiplier = multiplier * eta;
            end
            
            % do the necessary if the point is far enough
            if ( norm(X_k - Wnext ) > clickbacktol )
                expandsteps = expandsteps + 1;
                feascheckonlyFLAG =0;
                [~,nextxbar,nextx2bar,~] = simulate(y_seed,M_k,Wnext,feascheckonlyFLAG,OracleName);
                calls = calls + M_k;
                totlicalls = totlicalls + M_k;
                if (calls >= LIbudget)
                    budgetexceed_FLAG = 1;
                    return;
                end
                % update function value at candidate point; if we update global
                % list first, we will only be using partial crn
                %[~, nextxbar, ~] = updategloballist(Wnext,M_k,nextxbar,nextx2bar,iteration);
                Fnext = nextxbar;
                % update if Wnext is better
                if ( Fnext < F_k )
                    X_k = Wnext;
                    F_k = Fnext;
                    % update step size
                    l = step_factor * l;
                else
                    % stop the expansion phase
                    expandFLAG = 0;
                end
            else
                % stop the expansion phase because you've reached an infeasible pt
                expandFLAG = 0;
            end
        end
        
        % perform the click-back
        multiplier = sqrt(d);
        while ( norm(X_k - Wnext ) > clickbacktol)
            % identify Wnext
            %Wnext = round( 0.5 * ( X_k + Wnext ) );
            Wnext = round( X_k + ( multiplier * dhat_k) );
            % simulate at the next point if feasible
            feascheckonlyFLAG =0;
            [~,nextxbar,~,~] = simulate(y_seed,M_k,Wnext,feascheckonlyFLAG,OracleName);
            % do the necessary if the point is feasible
            if ( nextxbar ~= inf )
                calls = calls + M_k;
                totlicalls = totlicalls + M_k;
                if (calls >= LIbudget)
                    budgetexceed_FLAG = 1;
                    return;
                end
                % update function value at candidate point; if we update global
                % list first, we will only be using partial crn
                %[~, nextxbar, ~] = updategloballist(Wnext,M_k,nextxbar,nextx2bar,iteration);
                Fnext = nextxbar;
                % update if Wnext is better
                if ( Fnext < F_k )
                    % swap Wnext and X_k
                    Wnexttemp = X_k;
                    X_k = Wnext;
                    Wnext = Wnexttemp;
                    F_k = Fnext;
                end
            end
            multiplier = 0.9 * multiplier;
        end        
          
        % update number of line searches
        N_k = N_k + 1;
        % get out if you have done enough line searches or the line search
        % was very short
        if ( ( N_k >= b ) || ( N_k >=2 & expandsteps == 1 ) ) 
            timetoperformNEFLAG = 1;
            return;
        end

        % perform DA for another line search
        fromneFLAG = 0; Xtilde_k = X_k;
        if ( (N_k < b) && timetoperformNEFLAG == 0 )
            [daseed,y_seed, dhat_k, X_k, M_k, DAcalls,timetoperformNEFLAG] = DA(daseed,y_seed,X_k, Xtilde_k, F_k, M_k, lambda_k,iteration,LIbudget - calls,fromneFLAG,OracleName);
            calls = calls + DAcalls;
            totdacalls = totdacalls + DAcalls;
            % store the max sample size along nonlinear line search
            if M_k > maxSampleSize
                maxSampleSize = M_k;
            end
            if (calls >= LIbudget)
                budgetexceed_FLAG = 1;
                return;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_seed,Xtilde_k,betternbrfoundFLAG,calls] = R_A(y_seed,X_k,F_k,M_k,RAbudget,iteration)
    global budgetexceed_FLAG
    % initialize
    countercheck = 1;
    betternbrfoundFLAG = 0;
    calls = 0;
    Xtilde_k = X_k;
    % check if there is a better neighbor of X_k at sample size M_k
    d = size(X_k,1);
    [delnbdX_k] = find_delnbdN1(X_k);
    no_delnbdX_k = size(delnbdX_k, 2);

    while( ( countercheck <= no_delnbdX_k ) && ( betternbrfoundFLAG == 0 ) )
        nextpt = delnbdX_k(:,countercheck);
        feascheckonlyFLAG = 0;
        [~,nextxbar,nextx2bar,~] = simulate(y_seed,M_k,nextpt,feascheckonlyFLAG,OracleName);
        calls = calls + M_k;
        if (calls >= RAbudget)
            budgetexceed_FLAG = 1;
            return;
        end
        % update function value in the global list
        [~, nextxbar, ~] = updategloballist(nextpt,M_k,nextxbar,nextx2bar,iteration);
        if ( nextxbar < F_k )
            betternbrfoundFLAG = 1;
            Xtilde_k = nextpt;
        end
        countercheck = countercheck + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
