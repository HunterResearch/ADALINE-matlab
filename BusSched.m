%function [flag1, flag2, fn, FnVar, FnGrad, FnGradCov, constraint, ConstraintCov, ConstraintGrad, ConstraintGradCov, iseed]=BusSched(param, x, m, iseed, ri)
function [nextseed, y]=BusSched(nextseed, x)
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
% initialize 
param = [2 9 1 0 10 100];
% x=[10 20 30 40 50 60 70 80 90]
m=1;
iseed = nextseed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn=NaN;
FnVar=NaN;
FnGrad = NaN;
FnGradCov = NaN;
constraint = NaN;
ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;


id     = int32(param(2));
nseeds = int32(param(3));
lambda = param(5);
gamma  = param(6);

    %fprintf(1, 'id=%d, nseeds=%d, lambda=%0.12f, gamma=%0.12f\n', id, nseeds, lambda, gamma);
    
% model parameters feasibility check
flag1 = 0;
if lambda < 0 || gamma < 0 || id~=length(x) || nseeds~=length(iseed) 
    flag1 = 1;
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
        [iseed, u]= u16807d(iseed);
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
%FnVar=((sum2/m)-fn^2)/m;
%fprintf(1, 'ri=%d, m=%d, fn = %.12f, FnVar=%.12f\n',ri, m, fn, FnVar);
end



function [iseed,u16807d]=u16807d(iseed)
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
u16807d=0;
while (u16807d<=0 || u16807d>=1)
    iseed = mod (iseed * 16807,2147483647);
    u16807d = iseed / 2147483648;
end
end
