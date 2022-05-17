function optgappct = find_optgap(x,OracleName)
    if strcmp(OracleName,'BusScheduling') == 1
        sum2 = 0;
        lambda = 10;
        gamma = 100;
        sortedx = sort(x);
        d = size(x,1);
        for j = 1 : ( d + 1 )
            if ( ( j > 1 ) && ( j < (d + 1) ) )
                sum2 = sum2 + ( sortedx (j) - sortedx(j-1) )^2;
            elseif ( j == 1 )
                sum2 = sum2 + ( sortedx(j) - 0 )^2;
            else
                sum2 = sum2 + ( gamma - sortedx(j-1) )^2; 
            end
        end
        f_x = ( lambda / 2 ) *  sum2;
        opt_value = lambda * gamma^2 / ( 2 * ( d + 1 ) );
        optgappct = 100 * ( f_x / (opt_value) ) - 100;
    elseif strcmp(OracleName,'DiscreteQuadratic') == 1
        % this code computes the optimality gap at a given x for the discrete quadratic 
        % problem 
        optgappct = norm(x);
    elseif strcmp(OracleName,'DynamicNews') == 1
        % this code computes the optimality gap at a given x for the discrete quadratic 
        % problem; for reference, the optimal location is to put all your stock in the second location
        n = length(x);
        a = 10;
        b = 5 * n;
        c0 = 0.1;
        p = 0.5;
        xstarlast = floor(b- (c0/p) * ( b - a + 1 ));
        xstar = ([zeros(1,n-1),xstarlast])';
        optgappct = norm(x-xstar) / norm(xstar);
        %optgappct = 100 * max(abs(x-xstar)) / max(xstar);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

