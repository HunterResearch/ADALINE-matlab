function [nextseed, y, feasibilityFLAG] = OracleDynamicNews(nextseed, x, feascheckonlyFLAG)
%function [fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov] = DynamicNews(x, runlength, seed, other)
% x is the row vector, quantity to buy of each product
% runlength is the number of days of demand to simulate
% seed is the index of the substreams to use (integer >= 1)
% other is not used
% Returns Mean and Variance of Profit

%   *************************************************************
%   ***             Written by Danielle Lertola               ***
%   ***         dcl96@cornell.edu    June 27th, 2012          ***
%   ***              Edited by Bryan Chong                    ***
%   ***        bhc34@cornell.edu    October 15th, 2014        ***
%   ***             Edited by Raghu Pasupathy                 ***
%   ***        pasupath@purdue.edu  May 13, 2021              ***
%   *************************************************************
iseed = nextseed;
runlength = 1;
constraint = NaN;
ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;
FnGrad=NaN;
FnGradCov=NaN;
y=NaN;

% check the feasibility of x
feasibilityFLAG = 1;
if ( min(x) < 0 ) 
    feasibilityFLAG = 0;
    y = inf;
end

% if you are just checking feasibility, get out
    if ( feascheckonlyFLAG == 1 )
        return;
    end

if (max(x < 0)>0) | (runlength <= 0) | (runlength ~= round(runlength)) | (iseed <= 0) | (round(iseed) ~= iseed),
    %fprintf('All values in x should be >= 0, runlength should be positive integer, seed must be a positive integer\n');
    y = NaN;
    FnVar = NaN;
    FnGrad = NaN;
    FnGradCov = NaN;
else
    %%% Setting 1 %%%
    
    n=size(x,1);    % number of products
    % the number of customers is random
    [iseed, u]= mrg32k3a(iseed);
    a = 10;
    b = 5 * n;
    T = ceil( a + ( (b-a) * u ) ); % number of customers
    vnot=ones(n,1); % product constant
    for i=1:n
        vnot(i)=vnot(i)+i;
    end
    mu=100;
    c0 = 0.1;         % unit storage cost
    
    %%%%%%%%%%%%%%%%%
    
    
    % %%% Setting 2 %%%
    %
    % n=10; % number of products
    % T=30; % number of customers
    % vnot=ones(n,1)*5; % product constant
    % for i=1:n
    %    vnot(i)=vnot(i)+i;
    % end
    % mu=1;
    %
    % %%%%%%%%%%%%%%%%%
    
    % the i-th product's unit profit is (i-1)/n * 0.5, crowding the
    % products more and more as the number increases; the highest unit
    % profit is p = 0.5
    cost = (ones(n,1))';
    %sellPrice = (1:n);
    p = 0.5;
    sellPrice = cost + (0 : (n-1)) * (p / (n-1));
    
    % Generate a new stream for random numbers
    %OurStream = RandStream.create('mrg32k3a');
    
    % Set the substream to the "seed"
    %OurStream.Substream = nextseed;
    
    % Compute Gumbel RV's for Utility
    %OldStream = RandStream.setGlobalStream(OurStream);
    %Gumbel= evrnd(mu*-psi(1),mu,[n,runlength,T]);
    %RandStream.setGlobalStream(OldStream);
    
    Gumbel = zeros(n,runlength,T);
    for i = 1:n
        for run = 1 : runlength
            for j = 1:T
                [iseed, u]= mrg32k3a(iseed);
                Gumbel(i,run,j) = ( mu*-psi(1) ) - ( mu * log( log(1/u) ) );
            end
        end
    end
    
    % Determine Utility Function
    Utility=zeros(n,runlength,T);
    for i=1:n
        Utility(i,:,:)=vnot(i)+ Gumbel(i,:,:);
    end
    
    % Run Simulation
    initial=x*ones(1,runlength);
    inventory=initial;
    
    for j=1:T
        available=(inventory>0);
        decision=available.*Utility(:,:,j);
        [maxVal, index]=max(decision);
        itembought=maxVal>0;
        for k=1:runlength
            inventory(index(k),k)=inventory(index(k),k)-itembought(k);
        end
    end
    
    % Compute negative daily profit
    numSold =initial - inventory;
    unitProfit=sellPrice-cost;
    singleRepProfit=unitProfit*numSold;
    y = -mean(singleRepProfit) + c0 * sum(x);
    FnVar = var(singleRepProfit)/runlength;
    nextseed = iseed;
end