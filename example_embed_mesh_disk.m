%initialization
init
%load some data
load bump
%compute the embedding
embedder=embed_mesh(V,T,inds,true);
figure;
subplot(1,2,1);
title('original mesh');
embedder.visualize('dim',3);
campos([ -3.8165  -17.9552   11.5296]);
subplot(1,2,2);
title('embedding');
embedder.visualize('dim',2);
draw_disk
