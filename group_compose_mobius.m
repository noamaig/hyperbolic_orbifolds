function [ Ms ] = group_compose_mobius( M )
%given set of mobius trans, compute all their composisitions
a=[];b=[];c=[];d=[];
for i=1:length(M)
    a=[a;M{i}.a];
    b=[b;M{i}.b];
    c=[c;M{i}.c];
    d=[d;M{i}.d];
end
[I,J]=meshgrid(1:length(M),1:length(M));
I=I(:);J=J(:);
comp(:,1)=a(J).*a(I)+b(J).*c(I);
comp(:,2)=a(J).*b(I)+b(J).*d(I);
comp(:,3)=c(J).*a(I)+d(J).*c(I);
comp(:,4)=c(J).*b(I)+d(J).*d(I);
for i=1:3
    comp(:,i)=comp(:,i)./comp(:,4);
end
comp(:,4)=ones(size(comp,1),1);
tol=1e-8;
temp_comp=[real(comp) imag(comp)];
temp_comp = builtin('_mergesimpts',temp_comp,[tol tol tol tol tol tol tol tol],'first');
comp=temp_comp(:,1:4)+1i*temp_comp(:,5:8);
Ms=M;
for i=1:size(comp,1)
    Ms{end+1}=MobiusTrans(comp(i,1),comp(i,2),comp(i,3),comp(i,4));
end
end

