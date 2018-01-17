function [ M,P ] = hyperbolic_automorphism(  )
%the singularities
sings=[4 6 14];
assert(all(sings(1:end-1)<sings(2:end)));
%compute the angles from the singularities
angs=2*pi./sings;
%compute the arc lengths of the hyperbolic tri
h_lens=hyperbolic_lengths_from_angles(angs);
%only need the arcs that touch (0,0)
h_lens=h_lens(2:3);
%transform the lengths from hyperbolic to actual euclidean lengths
lens=hyperbolic_distance_to_euclidean_distance(h_lens);
%compute the positioning of the points 
p1=[lens(2) 0];
p3=[-lens(2) 0];
p2=[0 lens(1)];
p=[p1;p2;p3];
%transforming R^2 to complex numbers
p=p(:,1)+1i*p(:,2);
%these are the angles of the triangle from p1,p2,p3 - they are different
%than the angles used earlier because we constructed half a triangle (.e.g,
%a "white" triangle in the orbifold diagarms in wikipedia) and then glued
%it with a reflection when we chose p1,p2,p3 (i.e., a white+black tri in
%the diagrams).
real_angs=[angs(1)*2 angs(2) angs(3)*2];
%computing mobius trans which are automorphisms of the oribfold
%rotation around the p(2) 
M1=mobius_rotate_around_point(p(2),real_angs(3));
%rotation of pi around (0,0)
 M2=MobiusTrans(-1,0,0,1);
 %traversing the vertices of the polygon counterclockwise
 P=[p3;0 0;p1;p2]
 M={nan, M2,nan,M1};
 

end

