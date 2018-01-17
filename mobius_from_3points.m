function res=mobius_from_3points(s,t)
%%% compute the mobius trans taking 3 points to 3 points. s and t are each 3x2 array of 3 2d points

%f(z)=(az+b)/(cz+d)
%entails
%as_i+b=(cs+d)t_i
%s_i*a+1*b-s_it_i*c-t_id=0
s=s(:,1)+1i*s(:,2);
t=t(:,1)+1i*t(:,2);
M=[s ones(3,1) -s.*t -t; 1 0 0 0];
res=M\[0 0 0 1i+1]';
end