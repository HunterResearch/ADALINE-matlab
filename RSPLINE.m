%==========================================================================
%                        The R-SPLINE Algorithm 
%==========================================================================
% DATE
%        November 2014
%
% AUTHOR
%        Kalyani Nagaraj
%        kalyanin AT purdue DOT edu
%
% REFERENCE
%        H. Wang, B. Schmeiser, and R. Pasupathy. ACM TOMACS. 2013
%==========================================================================
%
% INPUT
%        x0 
%              Matrix (size = 'dim' X 'NumStartingSol') of 'NumStartingSol' 
%              initial solutions to the solver
%              For R-SPLINE, NumStartingSol=1
%              Each initial solution is of size dim X 1 
%        problem
%              Problem function name 
%        problemseed
%              Substream index (integer >=1)
%        solverseed
%              Input seed for R-SPLINE (integer between 1 and 2147483646)
%              See function u16807d for details.
%        budget
%              Vector of size NumSoln, where NumSoln is 
%              the number of solutions returned by the solver
%              for example, if budget = [500 1000] then NumSoln
%              is 2 and the solver returns best available solutions after 
%              every 500 calls to the oracle
%        logfilename       
%              
%
% OUTPUT
%        Ancalls
%              An array (size = 'NumSoln' X 1) of buget expended 
%        A 
%              An array (size = 'NumSoln' X 'dim') of solutions
%              returned by solver
%        Afn   
%              An array (size = 'NumSoln' X 1) of estimates of expected 
%              objective function value
%        AFnVar
%              An array of variances corresponding to  
%              the objective function at A
%              Equals NaN if solution is infeasible
%        AFnGrad
%              An array of gradient estimates at A; not reported
%        AFnGardCov
%              An array of gradient covariance matrices at A; not reported
%        Aconstraint
%              A vector of constraint function estimators; not applicable
%        AConstraintCov
%              An array of covariance matrices corresponding to the
%              constraint function at A; not applicable
%        AConstraintGrad
%              An array of constraint gradient estimators at A; not
%              applicable
%        AConstraintGradCov
%              An array of covariance matrices of constraint gradient  
%              estimators at A; not applicable
%
%==========================================================================

%% RSPLINE
function [alg_path] = RSPLINE(iseed, x0, budget,OracleName)

global pts_visited
global budgetexceed_FLAG
pts_visited=[];
 
% Solver parameters
mk=2;       % initial sample size (for k=1)
bk=10;      % mazximum number of SPLINE calls in the first retrospective 
            % iteration 
c1=1.1;     % growth rate of mk
c2=1;     % growth rate of bk

%logfname=strcat(logfilename,'.txt');
%logfid=fopen(logfname, 'w');

ncalls=0;	% tracks the total calls made to the oracle
x1=x0;
k=1;
alg_path = [[0,x0'],find_optgap(x0,OracleName)];

% use iseed(1) to create y_seed
seed = iseed(1,1); problemseed = zeros(6,1);
for j = 1 : 6
    [seed, u] = u16807d( seed );        
    problemseed(j,1) = floor(10^6 * u );
end


%Begin Retrospective Iterations
while ncalls<budget

    mk=ceil(mk*c1);
    bk=ceil(bk*c2);
    
    %fprintf(logfid, 'mk = %d, bk = %d\n', mk, bk);
    
    xk=x1;

	[splinencalls, x1, ~, ~] = SPLINE(xk, mk, bk, problemseed,OracleName);
    f_x = find_optgap(x1,OracleName);
    if splinencalls == -1     % initial solution is infeasible!
		return
    end
    ncalls = ncalls + splinencalls; %splinencalls is the number of oracle 
                                    %calls for each call of SPLINE
	alg_path = [alg_path; [ncalls, x1',f_x]];
    %fprintf('%4d %6d %10.3f %10.3f %10.3f %10d \n',k, mk, f_x, x1, ncalls);
    k=k+1;
    
    % advance the start seed if necessary to perform independent sampling 
    advanceby = 2*mk;
    j = 0;
    while ( j < advanceby ) 
        [problemseed, ~] = mrg32k3a(problemseed);
        j = j + 1;
    end
end
%fclose(logfid);

end

%% SPLINE
function [ncalls, xnew, xnewfn, xnewFnVar]=SPLINE(x0, mk, bk, iseed,OracleName)

ncalls=0;
xnew=x0;

[~,xnewfn,xnewFnSum2bar,~] = simulate(iseed,mk,xnew,0,OracleName);
xnewFnVar = ( xnewFnSum2bar - xnewfn^2 ) / mk;
ncalls = ncalls + 1;

if isnan(xnewfn) % infeasible solution
    ncalls=-1;
    return
end
x0fn=xnewfn;
x0FnVar=xnewFnVar;

%for i=1:bk
while ( ncalls <= bk )
		
	[SPLIncalls, xold, xoldfn, xoldFnVar] = SPLI(xnew, xnewfn, xnewFnVar, mk, iseed,OracleName);
		
	[NEncalls, xnew, xnewfn, xnewFnVar] = NE(xold, xoldfn, xoldFnVar, mk, iseed,OracleName);

    ncalls = ncalls + SPLIncalls + NEncalls;
        
    if  xoldfn==xnewfn
		%fprintf('\n\t\tSPLINE ended at bk=%d since NE and SPLI returned the same solution\n\n', i);
		break
    end
end    
	
% starting solution is better than the new solution by a very small margin
if x0fn + 0.00005 <= xnewfn
	xnew=x0;
	xnewfn=x0fn;
    xnewFnVar=x0FnVar;
end
end

%% NE
function [ncalls, xnew, xnewfn, xnewFnVar] = NE(x, fn, FnVar, mk, iseed,OracleName)

id = length(x); 
xold=x;
xnew=x;
xnewfn=fn;
xnewFnVar=FnVar;

ncalls=0;
y2=fn;
ixquad=zeros(id,1);

for i=1:id
	count=1;

    xold(i)=xold(i)+1;
    [~,xoldfn, xoldFnSum2bar,~]=simulate(iseed, mk, xold,0,OracleName);
    xoldFnVar = ( xoldFnSum2bar - xoldfn^2 ) / mk;

    if ~isnan(xoldfn)
		ncalls = ncalls + mk;
		y1=xoldfn; 
		count=count+1;
        if xoldfn + 0.00005 < xnewfn
			xnew=xold;
			xnewfn=xoldfn;
			xnewFnVar=xoldFnVar;
        end
    end

	xold(i)=xold(i)-2;
    [~,xoldfn, xoldFnSum2bar,~]=simulate(iseed, mk, xold,0,OracleName);
    xoldFnVar = ( xoldFnSum2bar - xoldfn^2 ) / mk;
    if ~isnan(xoldfn)
		ncalls = ncalls + mk;
		y3=xoldfn; 
		count=count+1;
        if xoldfn + 0.00005 < xnewfn
			xnew=xold;
			xnewfn=xoldfn;
            xnewFnVar=xoldFnVar;
        end 
    end
	xold(i)=xold(i)+1;
	xqnew=xold(i);

    %quadratic search
    if count==3 
		a = (y1+y3)/2.0 - y2;
		b = (y1-y3)/2.0;
        if a-0.00005 > 0
			xqnew = int32(xold(i) - (b / (a + a)));
        end
    end
    if  abs(xqnew) < 2147483646.0 %2^31-2
		ixquad(i) = xqnew;
    end
end
	
[~,ixquadfn, ixquadFnSum2bar,~]=simulate(iseed, mk, ixquad,0,OracleName);
ixquadFnVar = ( ixquadFnSum2bar - ixquadfn^2 ) / mk;

%if ~isnan(ixquadfn)
%	ncalls = ncalls + mk; 
%    if ixquadfn + 0.00005 < xnewfn
%		xnew=ixquad;
%        xnewfn=ixquadfn;
%        xnewFnVar=ixquadFnVar;
%    end
%end

end


%% SPLI
function [ncalls, xbest, xbestfn, xbestFnVar] = SPLI(x, fn, FnVar, mk, iseed,OracleName)

id = length(x); 
imax=100;
jmax=1;

xbest=x;
xbestfn=fn;
xbestFnVar=FnVar;

ncalls=0;
s0=2.0;
c=2.0;
			
for j=0:jmax
    %PRINT TO LOGFILE

    [x1, ~]=PERTURB(xbest, iseed);

   [~, gamma, npoints, plixbest, plixbestfn, plixbestFnVar] = PLI(x1, mk, iseed,OracleName);

	ncalls = ncalls + npoints*mk;	

	%regardless of whether npoints=id+1 or not, update current best
    if  plixbestfn + 0.00005 < xbestfn && npoints>0
        xbest=plixbest;
        xbestfn=plixbestfn;
        xbestFnVar=plixbestFnVar;
    end
		
    if npoints < id+1
        return 
    end
			
	glength=norm(gamma);
		
    if glength + 0.00005 <= 0
		return
    end
    x0=xbest;
    gamma=gamma/glength;

    for i=0:imax
        s = s0 * c^i;
        ix1=zeros(id,1);
        for k=1:id
			ix1(k)=floor(x0(k)-s*gamma(k)+0.5);
        end
        [~,ix1fn, ix1FnSum2bar,~]=simulate(iseed, mk, ix1,0,OracleName);
        ix1FnVar = ( ix1FnSum2bar - ix1fn^2 ) / mk;

        if isnan(ix1fn) % if ix1 is infeasible
            return	
        end
		ncalls = ncalls + mk;
        if ix1fn >= xbestfn + 0.00005 && i <= 2
            return
        end
        if ix1fn >= xbestfn + 0.00005
            break
        end
        
		xbest=ix1;
		xbestfn = ix1fn;
		xbestFnVar=ix1FnVar;
    end    
end
end


%% PLI
function [fbar, gamma, npoints, plixbest, plixbestfn, plixbestFnVar] = PLI(x, mk, iseed,OracleName)

id = length(x); npoints=0;
strange=3.145962987654;
gamma=zeros(id,1);
x0=floor(x);
%z=zeros(id,1);
z=x-x0;
z=[1;z;0];

[~, p]=sort(z, 'descend');
w=zeros(id+1,1);
for i=1:id+1		
	w(i,1)=z(p(i))-z(p(i+1));
end	
wsum=0;
fbar=0;

[~,x0fn, x0FnSum2bar,feasibleFLAG]=simulate(iseed, mk, x0,0,OracleName);
x0FnVar = ( x0FnSum2bar - x0fn^2 ) / mk;
if ( ~isnan(x0fn) && feasibleFLAG == 1 )
	npoints=npoints+1;
	wsum = wsum + w(1);
	fbar = fbar + w(1)*x0fn;
	ghatold = x0fn;
    plixbest=x0;
    plixbestfn = x0fn;
    plixbestFnVar=x0FnVar;	
    
else
    ghatold = 0;
	plixbestfn = strange;
end

for i=2:id+1
    x0(p(i)-1)=x0(p(i)-1)+1;
    
    %call oracle at the other id points that form the simplex
    [~,x0fn, x0FnSum2bar,~]=simulate(iseed, mk, x0,0,OracleName);
    x0FnVar = ( x0FnSum2bar - x0fn^2 ) / mk;
    if ~isnan(x0fn)
		npoints=npoints+1;
		wsum = wsum + w(i);
		fbar = fbar + w(i)*x0fn;
		gamma(p(i)-1) = x0fn - ghatold;
		ghatold = x0fn;
		
        
        if plixbestfn == strange || x0fn + 0.00005 < plixbestfn
            plixbest=x0;
            plixbestfn = x0fn;
            plixbestFnVar=x0FnVar;	
        end
    end
end

if wsum > 0.00005
    fbar = fbar/wsum;	
end

end


%% PERTURB
function [xpert, sseed] = PERTURB(x, sseed)

id=length(x);
xpert=zeros(id,1);

for i=1:id
	[sseed, u] = mrg32k3a(sseed);
    xpert(i) = x(i) + .3*(u - 0.5);
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
