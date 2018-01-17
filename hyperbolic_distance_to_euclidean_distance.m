function [ lens ] = hyperbolic_distance_to_euclidean_distance( h_lens)
%transform distances (from 0,0 to some) in poincare to euclidean distances from
%(0,0) 

%H=ln((1+E)/(1-E))
%e^H=(1+E)/(1-E)
%e^H(1-E)=1+E
%E(1+e^H)=e^H-1
%E=(e^H-1)/(1+e^H)
lens=(exp(h_lens)-1)./(1+exp(h_lens));
end

