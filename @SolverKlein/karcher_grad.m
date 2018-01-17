function [ga,gb,gc,gd,hd]=karcher_grad(obj,a,b,c,d,w)
            %d(u,v) in hyperboloid is archosh(U./|U|.*V./|V|), where U.*V=u1.*v1-u2.*v2-u3.*v3
            %this is projected on klein onto by setting timlike component, u1 v1 to 1
            %so we get 
           
            uv=1- a.*c-b.*d;
            nu=sqrt(1-a.^2-b.^2);
            nv=sqrt(1-c.^2-d.^2);
            hd=acosh(uv./(nu.*nv));
            
 ga=  -(2.*acosh(-(a.*c + b.*d - 1)./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(1./2)).*(c./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(1./2) - ((2.*a.*c.^2 + 2.*a.*d.^2 - 2.*a).*(a.*c + b.*d - 1))./(2.*(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(3./2))))./((a.*c + b.*d - 1).^2./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1) - 1).^(1./2);
 gb=-(2.*acosh(-(a.*c + b.*d - 1)./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(1./2)).*(d./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(1./2) - ((2.*b.*c.^2 + 2.*b.*d.^2 - 2.*b).*(a.*c + b.*d - 1))./(2.*(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(3./2))))./((a.*c + b.*d - 1).^2./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1) - 1).^(1./2);
 gc=            -(2.*acosh(-(a.*c + b.*d - 1)./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(1./2)).*(a./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(1./2) - (c.*(a.*c + b.*d - 1).*(a.^2 + b.^2 - 1))./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(3./2)))./((a.*c + b.*d - 1).^2./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1) - 1).^(1./2);
 gd=            -(2.*acosh(-(a.*c + b.*d - 1)./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(1./2)).*(b./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(1./2) - (d.*(a.*c + b.*d - 1).*(a.^2 + b.^2 - 1))./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1).^(3./2)))./((a.*c + b.*d - 1).^2./(c.^2.*(a.^2 + b.^2 - 1) + d.^2.*(a.^2 + b.^2 - 1) - a.^2 - b.^2 + 1) - 1).^(1./2);
%next few lines are to handle when dividing by very small number when two
%vertices are very close and cause inf's and nan's
badinds=isnan(ga)|isinf(ga)|isnan(gb)|isinf(gb)|isnan(gc)|isinf(gc)|isnan(gd)|isinf(gd);
ga(badinds)=0;
gb(badinds)=0;
gc(badinds)=0;
gd(badinds)=0;
       ga=w.*ga;     
       gb=w.*gb;
       gc=w.*gc;
       gd=w.*gd;
            
            
        end
        