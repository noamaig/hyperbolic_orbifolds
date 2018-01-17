function computeDistortion(obj,force)

if nargin<2
    force=false;
end
if ~force && ~isempty(obj.dets)
    disp('using precomputed dist');
    return;
end
orgV=obj.flat_V;
obj.flat_V=obj.toKlein();
if isempty(obj.As)||force
    obj.computeAs();
end
fprintf('==== Computing distortion etc. ===\n');
tid=tic;
As=obj.As;


a = squeeze(As(1,1,:))';
b = squeeze(As(1,2,:))';
c = squeeze(As(2,1,:))';
d = squeeze(As(2,2,:))';%the entries of A

obj.dets=(a.*d-b.*c)';
obj.frobenius=sqrt(a.^2+b.^2+c.^2+d.^2)';


alpha = [a+d;b-c]/2; %2XM
beta = [a-d;b+c]/2;
alpha=sqrt(sum(alpha.^2));
beta=sqrt(sum(beta.^2));
obj.smax=(alpha+beta)';
obj.smin=abs((alpha-beta)');

obj.flipped=obj.dets<=0;
if any(obj.flipped)
    warning('there are flipped tris');
end
fprintf('done computing, %f seconds\n',toc(tid));
if any(obj.flipped)
    fprintf('!!!! FLIPPED: %d\n',nnz(obj.flipped));
end
obj.flat_V=orgV;
end
