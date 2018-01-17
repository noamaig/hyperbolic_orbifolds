function [grada2,gradb2,gradc2,gradd2,hd2]=karcher_grad(obj,a,b,c,d,w)
            
            
            f=2*(a.^2+b.^2+c.^2+d.^2-2*a.*c-2*b.*d);
            g=1-a.^2-b.^2-c.^2-d.^2+a.^2.*c.^2+a.^2.*d.^2+b.^2.*c.^2+b.^2.*d.^2;
            delta2=f./g;
            hd2=real((acosh(1+delta2)));
            f_a=4*a-4*c;
            f_b=4*b-4*d;
            f_c=4*c-4*a;
            f_d=4*d-4*b;
            g_a=2*a.*(c.^2+d.^2-1);
            g_b=2*b.*(c.^2+d.^2-1);
            g_c=2*c.*(a.^2+b.^2-1);
            g_d=2*d.*(a.^2+b.^2-1);
            
            grad_delta_a=(g.*f_a-f.*g_a)./(g.^2);
            grad_delta_b=(g.*f_b-f.*g_b)./(g.^2);
            grad_delta_c=(g.*f_c-f.*g_c)./(g.^2);
            grad_delta_d=(g.*f_d-f.*g_d)./(g.^2);
            
            fac=(sqrt(delta2.^2+2*delta2));
            grada2=2*w.*hd2.*grad_delta_a./fac;
            gradb2=2*w.*hd2.*grad_delta_b./fac;
            gradc2=2*w.*hd2.*grad_delta_c./fac;
            gradd2=2*w.*hd2.*grad_delta_d./fac;
            grada2(hd2==0)=0;
            gradb2(hd2==0)=0;
            gradc2(hd2==0)=0;
            gradd2(hd2==0)=0;
            
            %         ab=(1-a.^2-b.^2).^2;
            %         cd=(1-c.^2-d.^2).^2;
            %
            %         grada2=grada2.*ab;
            %         gradb2=gradb2.*ab;
            %
            %         gradc2=gradc2.*cd;
            %         gradd2=gradd2.*cd;
            
            %         hd2=real(hd2);
            
            
            
            
            
        end
        