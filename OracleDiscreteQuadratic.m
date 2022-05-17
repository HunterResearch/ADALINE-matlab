function [nextseed,y, feasibilityFLAG] = OracleDiscreteQuadratic(nextseed,x,feascheckonlyFLAG)
    y = nan;
    kappa = 20;
    % the problem is boundconstrained
    feasibilityFLAG = 1;
    %if max( abs(x) > 15 )
    %    feasibilityFLAG = 0;
    %end
    % if you are just checking feasibility, get out
    if ( feascheckonlyFLAG == 1 )
        return;
    end

    d = size(x,1);
    % simple quadratic function plus noise
    %rng(nextseed(1));
    % obtain the state of the random number generator
    s = rng;
    % take it to a specific stat for construction  of f_x
    rng(788495);
    % generate the eigen values randomly between 1 and kappa + 1
    rc = 1 + ( kappa * rand(d,1));
    B = sprandsym(d,0.5,rc);
    %B = eye(size(x,1));
    if ( round(x) == x )
        f_x = x' * B * x;
    else
        f_x = inf;
    end
    
    % restore the state of the generator
    rng(s);
    [nextseed,u1] = mrg32k3a(nextseed);
    [nextseed,u2] = mrg32k3a(nextseed);
    sigma2 = 5;
    sigma1 = 5;
    e2 = norminv(u1,0,sigma2);
    e1 = 0;
    
    %if x~= 0
    %    if norm(x) > sqrt(d)
    c_1 = norm(B*x);
    e1 = c_1 * sigma1 * norminv(u2,0,1);
    %    else
    %        c_1 = 0;
    %        e1 = 0;
    %    end
    %end
    y = f_x + e1 + e2;
    if isnan(y)
        fprintf("Hello")
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
