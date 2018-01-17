function [ K_V ] = poincare_to_klein( P_V )
%convert points from the Poincare model to the Klein model.
n=1+sum(P_V.^2,2);
K_V=2*[P_V(:,1)./n P_V(:,2)./n];
end


