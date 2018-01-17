function [ M ] = disc_mobius_translation( p,theta )
%compute mobius transformation that maps the unit disc into itself,
%translates the(complex) point p to (0,0) and possibly rotates by R if given
if nargin<2
    theta=0;
end
assert(length(p)<=2);
if length(p)==2
    p=p(1)+1i*p(2);
end
R=exp(1i*theta);
a=1;
b=-p;
c=conj(b);
d=1;

M=MobiusTrans(R.*a, R.*b, c, d);

end

