function [ h_lens ] = hyperbolic_lengths_from_angles( angs )
%given 3 angles, compute arc lengths of hyperbolic tri from 2nd hyperbolic cosine law

alpha=angs;
beta=angs([2:end 1]);
gamma=angs([3 1:2]);

%cos(alpha)=-cos(beta)cos(gamma)+sin(beta)sin(gamma)cosh(alpha/k)
%cos(alpha)+cos(beta)cos(gamma)=sin(beta)sin(gamma)cosh(alpha/k)
%(cos(alpha)+cos(beta)cos(gamma))/sin(beta)sin(gamma)=cosh(alpha/k)
h_lens=acosh((cos(alpha)+cos(beta).*cos(gamma))./(sin(beta).*sin(gamma)));
end

