%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optgappct = FindGapBusWaitingTime(x)
% this code computes the optimality gap at a given x for the d-bus scheduling 
% problem with lambda = 10, gamma = 100
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
    opt_value = lambda * gamma^2 / ( 2 );
    optgappct = 100 * ( f_x / (opt_value)) - 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%