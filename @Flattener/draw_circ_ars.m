function hsteps=draw_circ_ars(obj)
E=edges(obj.flat_T);
I=E(:,1);
J=E(:,2);
pI=obj.flat_V(I,1)+1i*obj.flat_V(I,2);
pJ=obj.flat_V(J,1)+1i*obj.flat_V(J,2);


a=1;
b=-pI;
c=conj(b);
d=1;

mP=(a.*pJ+b)./(c.*pJ+d);
steps=mP*linspace(0,1,20);

ai=d;
bi=-b;
ci=-c;
di=a;

nom=bsxfun(@times,steps,ai);
nom=bsxfun(@plus,nom,bi);
denom=bsxfun(@times,steps,ci);
denom=bsxfun(@plus,denom,di);
hsteps=nom./denom;
 
 hold on
 for i=1:size(hsteps,1)
 line(real(hsteps(i,:)),imag(hsteps(i,:)));
 pause(0.001)
hold off

end