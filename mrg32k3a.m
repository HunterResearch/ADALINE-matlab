function [seed,u] = mrg32k3a(seed)
%
%..........................................................................
%     Raghu Pasupathy     June 2012.                                   
%     
%   This is an implementation of Pierre L'Ecuyer's Random Number Generator,
%   MRG32K3A. ("Good Parameters and Implementations for Combined Multiple 
%   Recursive Random Number Generators", Operations Research 47, pp.
%   159-164, 1999.)
%..........................................................................
%     modified by Kyle Cooper     July 2018.
%     Modified to match L'Ecuyer's implementation.
%..........................................................................

s1 = 1403580;
t1 = 810728;
s2 = 527612;
t2 = 1370589;
m1 = 4294967087;  % = 2^32 - 209
m2 = 4294944443;  % = 2^32 - 22853
m3 = 4294967088;  % = 2^32 - 208

%  Pasupathy version
%p1 = mod( ( s1 * seed(1) ) - ( t1 * seed(2) ), m1);
%p2 = mod( ( s2 * seed(4) ) - ( t2 * seed(5) ), m2);

% to match L'Ecuyer -- change seed indices and mod(x,a) to mod(mod(x, a), a)
p1 = mod(mod(s1 * seed(2) - t1 * seed(1), m1), m1);
p2 = mod(mod(s2 * seed(6) - t2 * seed(4), m2), m1);

z = mod( ( p1 - p2 ), m1 );

if ( z > 0 )
    u = z / m3;
elseif ( z == 0 )
    u = m1 / m3;
end

seed(1) = seed(2); seed(2) = seed(3); seed(3) = p1;
seed(4) = seed(5); seed(5) = seed(6); seed(6) = p2;

%--------------------------------------------------------------------------