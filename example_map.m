%% loading some data

load david_50k
embedder_david=Flattener(V,T,inds_for_map);
load maxp
embedder_max=Flattener(V,T,inds_for_map);
%% The main mapping part

%generating a cell-array of all the meshes we want to embed - this can be
%of any length > 1, and will define a simultaneous mapping between all
%meshes in the collection.
flatteners={embedder_david,embedder_max};
mapper=Mapper(flatteners);

%next line embeds all the meshes (if they haven't been already embedded),
%and computes the induced mapping of mesh 1 to all other meshes. So in this
%example computes the mapping from 1 to 2.
mapper.computeMap();
%... and this is how, given that the mapper already embedded all given meshes, 
% you compute the mapping from other pairs in the collection,
% in this case computing mapping from 2 to 1
mapper.lift(2,1);
% The result of the computation is in mapper.map.barCoords which is a
% cellarray of size NxN where N is the number of meshes in the collection
% (in our case =2). The (i,j) entry is the map i->j, or empty if it wasn't 
% computed, and then you need to call lift(i,j) to compute it. The map is 
% represented by a matrix of size |V_j| x |V_i| where V_k are the vertices 
% of the k'th mesh. 
% IMPORTANT DISCLAIMER: This is NOT the exact map defined by the embeddings.
% The exact mapping is a homeomorphism between the two surfacrs, in general
% mapping vertices of mesh 1 onto faces of mesh 2, and hence the
% barycentric mapping is just a rough estimation of the map for
% visualization, basically sampling the mapping only on the vertices.
% 

%% Visualization of the map as a movie with linear interpolation between the 2 models.

% next line is mapping all of max's vertices onto davis's mesh.
% note: this is confusing, the "geometric" map is max->david, but the 
% linear map of the barycentric system recieves as input the position of
% david's vertices and position's max's vertices according to it, so the
% input\output of geometric map vs. linear map is reversed. 
V_on_david=mapper.map.barCoords{1,2}*embedder_david.M_orig.V';


% Here I am just rotating the two point clouds to be kind of aligned so the
% linear interpolation will look decent. 
V_on_max=embedder_max.M_orig.V';
%center
V_on_david=bsxfun(@minus,V_on_david,mean(V_on_david));
V_on_max=bsxfun(@minus,V_on_max,mean(V_on_max));



%find R minimizing |A*R-B|=|AA|+|BB|-tr(R^T*A^T*B) -> maximize tr(R^TA^TB)
A=V_on_david'*V_on_max;
[U,~,V]=svd(A);
R=U*V';
V_on_david=V_on_david*R;
% min_a |a*A-B| -> min_a a^2*|A^2|-a*tr(A'*B)
% 2a*|A^2|=tr(A'B);
scale=trace(V_on_david'*V_on_max)/sum(V_on_david(:).^2);
V_on_david=V_on_david*scale;
%the movie with linear interpolation. Spoilers: it ends with david. 
figure(1);
clf;
set(gcf,'name','linear interpolation between source and target');
T=embedder_max.M_orig.F';
inc=0.05;
t=0;
hP=patch('faces',T,'vertices',V,'facecolor','white','edgealpha',0.3);
axis equal
light
campos([101.7697   57.8614  -97.7677]);
camup([0 1 0]);
while(true)
   t=t+inc;
   if t>=1
       t=1;
       inc=-inc;
   end
   if t<=0
       t=0;
       inc=-inc;
   end
   V=t*V_on_max+(1-t)*V_on_david;
   hP.Vertices = V;
   drawnow;
   %pause(0.05);
end