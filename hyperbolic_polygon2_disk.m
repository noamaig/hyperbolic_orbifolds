function [ realP,reflectionPoints] = hyperbolic_polygon2_disk(  k)

%we first treat the triangle T, which has one vertex o located at (0,0)
%with angle 2pi/2k, a vertex v with angle 2pi/2k and a vertex r with 90
%degrees. taking this triangle and reflecting along the edge (o,e) gets us
%the basic triangle T' we will use to construct the prefect polygon. T' has
%vertices o,v,v' with angles 2pi/k 2pi/2k 2pi/2k appropriately. After this
%gluing the vertex r has an angle of 180.

%the angles of T, according to o,v,r
angs=[2*pi/(k*2) pi/4 pi/2];
%compute the hyp arc lengths of T
h_lens=hyperbolic_lengths_from_angles(angs);
%only need the arcs that touch o
lens=hyperbolic_distance_to_euclidean_distance(h_lens);
%place the two points
v=[0 lens(3)]';
r=[0 lens(2)]';
theta=2*pi/(2*k);
R=[cos(theta) -sin(theta);
    sin(theta) cos(theta)];
r=R*r;
%transform the lengths from hyperbolic to actual euclidean lengths
h_lens=h_lens(1);
lens=hyperbolic_distance_to_euclidean_distance(h_lens);
%compute the positioning of the points

P=[];
M={};
for i=1:k
    theta=2*pi*i/k;
    R=[cos(theta) -sin(theta);
        sin(theta) cos(theta)];
    cur_v=R*v;
    cur_r=R*r;
    P=[P cur_v cur_r];
    M1=mobius_rotate_around_point(cur_r(1)+1i*cur_r(2),pi);
    %rotation of pi around (0,0)
    M{end+1}=nan;
    M{end+1}=M1;
    
end
P=P(1,:)+1i*P(2,:);
P=P';

% P is the original polygon, Q is the duplicated
P=[real(P) imag(P)];
realP=P(1:2:end,:);
reflectionPoints=P(2:2:end,:);
end

