function P_V= klein_to_poincare( K_V )
%convert points in Klein model to Poincare model
n=1+sqrt(1-sum(K_V.^2,2));
P_V=K_V;
for i=1:2
P_V(:,i)=P_V(:,i)./n;
end


end

