function [Ti,B,BC]=lift(obj,other)
tid=tic;
%transfer all vertices to klein
myV=obj.toKlein();

%the automorphisms
if ~obj.isdisc
    M=automorphisms(obj.M,obj.P);
    M=group_compose_mobius( M );
else
    M={};
end

%poincare vertices of other mesh
orgOtherV=other.flat_V;
%generate my triangulation in klein
TR=triangulation(obj.flat_T,myV);
%add the identity transformation to the automorphisms
M=[{MobiusTrans(1,0,0,1)} M];
%vertices of other mesh still not found
unmatched=1:length(orgOtherV);
%bar coords
B=zeros(length(unmatched),3)-1;
%triangle in source matching target
Ti=zeros(length(unmatched),1)-1;
%we already know matching of cones
unmatched=setdiff(unmatched,other.cut_cone_inds);


if obj.isdisc
    otherV=other.toKlein();
    for j=1:length(obj.reflection_paths)
        myp=obj.reflection_paths{j};
        otherp=other.reflection_paths{j};
        
        n=normr(-myV(myp(1),:)+myV(myp(end),:));
        myVals=myV(myp,:)*n';
        otherVals=otherV(otherp,:)*n';
        for k=2:length(otherVals)-1
            v=otherVals(k);
            ei=find(v>=myVals(1:end-1)&v<myVals(2:end));
            ei=ei(1);
            e=[myp(ei) myp(ei+1)];
            ti=find(sum(ismember(obj.flat_T,e),2)==2);
            assert(length(ti)==1);
            otherind=otherp(k);
            Ti(otherind,:)=ti;
            B(otherind,obj.flat_T(ti,:)==e(1))=(v-myVals(ei))/(myVals(ei+1)-myVals(ei));
            B(otherind,obj.flat_T(ti,:)==e(2))=(v-myVals(ei+1))/(myVals(ei+1)-myVals(ei));
        end
        unmatched=setdiff(unmatched,otherp);
        %                 otherVKlein(p,:)=otherVKlein(p,:)+kron(other.reflections_shifts(i)-otherVKlein(p,:)*n',n);
    end
end

%now go over all automorphisms and find matches
for i=1:length(M)
    %         fprintf('Lifting: in automorphism #%d/%d, %d/%d left to lift\n',i,length(M),length(unmatched),length(orgOtherV));
    %apply automorphism to vertices of other mesh
    otherVPoincare=M{i}.map(orgOtherV(unmatched,:));
    %transfer the vertices to klein
    otherVKlein=poincare_to_klein(otherVPoincare);
    
    %find matches
    [ti,b] = TR.pointLocation(otherVKlein);
    found=~isnan(ti);
    if ~any(found)
        continue;
    end
    
    matched=unmatched(found);
    B(matched,:)=b(found,:);
    Ti(matched)=ti(found);
    unmatched=unmatched(~found);
    
    
    %now correct
    V=M{i}.map(orgOtherV);
    for j=1:length(matched)
        
        ind=matched(j);
        
        MM=disc_mobius_translation(V(ind,:));
        t=obj.flat_T(Ti(ind),:);
        t=obj.flat_V(t,:);
        t=MM.map(t);
        %             tr=triangulation([1 2 3],t);
        A=[t';1 1 1];
        b=A\[0;0;1];
        mm=min(b);
        if mm<-1e-8
            error('negative bar coords: %f %f %f',b);
        elseif mm<0
            b=b-mm;
            b=b/sum(b);
        end
        %             [ti,b] = tr.pointLocation([0 0]);
        B(ind,:)=b;
    end
    
    
    
    if isempty(unmatched)
        break;
    end
end

if ~isempty(unmatched)
    
    fprintf('not found inds: %d\n',unmatched);
    error('During lifting, couldn''t find image of some vertices : %d out of %d',length(unmatched),length(other.flat_V));
end
%completing bar coords for the cones
[~,tinds]=ismember(obj.cut_cone_inds,obj.flat_T);
% ii - the first tri in obj that has the cone
% jj - the position of the cone in the tri
[ii,jj]=ind2sub(size(obj.flat_T),tinds);
%matching triangles to the cones of other
Ti(other.cut_cone_inds)=ii;
%matching bar coords to the cones of other
B(other.cut_cone_inds,:)=sparse(1:length(jj),jj,ones(size(jj)),length(jj),3);
I=repmat((1:length(other.flat_V)),1,3);
I=I(:);
J=obj.flat_T(Ti,:);
J=J(:);
V=B(:);
BC=sparse(I,J,V,length(other.flat_V),length(obj.flat_V));
BC(sub2ind(size(BC),other.cut_cone_inds,obj.cut_cone_inds))=1;
obj.times.lift=toc(tid);
end