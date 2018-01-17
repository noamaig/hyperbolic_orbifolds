function [ M ] = mobius_rotate_around_point( p,theta )
%compute the mobius transformation that represents a rotation of angle
%theta around the point p (in the poincare model)
 M_T = disc_mobius_translation( p);

 a=exp(1i*theta);
 b=0;
 c=0;
 d=1;
 M_R=MobiusTrans(a,b,c,d);
 M_I=M_T.inverse();
 M=M_I.compose(M_R.compose(M_T));
end

