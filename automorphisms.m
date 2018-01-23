function  Ms = automorphisms( Morig,cones ,iterates)
%generate the automorphisms of the group associated with the orbifold
if nargin<3
    iterates=1;
end

if size(cones,2)==2
    c=cones(:,1)+1i*cones(:,2);
else
    c=cones;
end
for i=1:length(Morig)
    Morig{end+1}=Morig{i}.inverse();
end
Morig{end+1}=MobiusTrans(1,0,0,1);
for i=1:length(Morig)
    A(i,1:4)=[Morig{i}.a Morig{i}.b Morig{i}.c Morig{i}.d];
end
Ms=Morig;
for i=1:length(Morig)
    Mcur=Morig{i}.inverse();
    d=sum(abs(bsxfun(@minus,A,[Mcur.a Mcur.b Mcur.c Mcur.d])),2);
            if min(d)<1e-8
                continue;
            end
            Ms{end+1}+Mcur;
end
Ms=helper( Ms,c,A ,iterates);
end
function Ms=helper( Morig,c,A ,iterates)
if iterates==0
    Ms=Morig;
    return;
end
Ms=Morig;
for iter=1:1000
    curLen=length(Ms);
    for i=1:length(Morig)
        for j=1:length(Ms)
            Mcur=Morig{i}.compose(Ms{j});
            d=sum(abs(bsxfun(@minus,A,[Mcur.a Mcur.b Mcur.c Mcur.d])),2);
            if min(d)<1e-8
                continue;
            end
            d=abs(c-Mcur.map(c));
            if min(d)>1e-8
                continue;
            end
            Ms{end+1}=Mcur;
            A(end+1,1:4)=[Mcur.a Mcur.b Mcur.c Mcur.d];
        end
    end
    
    if curLen==length(Ms)
        break;
    end
end
newCones=c;
for i=1:length(Ms)
    newCones=[newCones;Ms{i}.map(c)];
end
% temp_comp = builtin('_mergesimpts',newcones,[tol tol tol tol tol tol tol tol],'first');

Ms=helper( Ms,newCones,A ,iterates-1);
end

