init
load david_50k
embedder=embed_mesh(V,T,inds,false);
figure;
subplot(1,2,1);
title('original mesh');
embedder.visualize('dim',3);campos([ -6.1034    9.7864    1.9207]);
subplot(1,2,2);
title('embedding');
embedder.visualize('dim',2);
draw_disk
