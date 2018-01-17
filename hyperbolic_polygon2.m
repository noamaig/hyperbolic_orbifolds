function [ FullP,M,matching_vertices,distance_matrix] = hyperbolic_polygon2(  k,symmetric)

%we first treat the triangle T, which has one vertex o located at (0,0)
%with angle 2pi/2k, a vertex v with angle 2pi/2k and a vertex r with 90
%degrees. taking this triangle and reflecting along the edge (o,e) gets us
%the basic triangle T' we will use to construct the prefect polygon. T' has
%vertices o,v,v' with angles 2pi/k 2pi/2k 2pi/2k appropriately. After this
%gluing the vertex r has an angle of 180.

if nargin<2
    symmetric=true;
end
if ~symmetric
    warning('using non-symmetric embedding');
end
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
% h_lens=h_lens(1);
% lens=hyperbolic_distance_to_euclidea0n_distance(h_lens);
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

Mr=mobius_rotate_around_point(P(2),pi);
Q=Mr.map(P);
% P=[P;Q(end:-1:1,:)];

% P is the original polygon, Q is the duplicated
P=[real(P) imag(P)];
Q=[real(Q) imag(Q)];

%setting Q so it is the exact revese of P and Q(i) matches P(i) in the
%orbifold
Q=[Q(4:end,:);Q(1:3,:)];
Q=Q(end:-1:1,:);

%set the first vertex to be first
P=[P(3:end,:);P(1:2,:)];
Q=[Q(3:end,:);Q(1:2,:)];


%last vertex is not relevant, remove it
P=P(1:end-1,:);
Q=Q(1:end-1,:);
s=P(1,:);
e=P(end,:);
av=(s+e)/2;
%reversing Q
FullP=[P;Q(end-1:-1:2,:)];
FullP=FullP(1:2:end,:);
%removing points on Q which are also on P

% P=P(1:2:end,:);
% Q=Q(2:2:end,:);
% P=P(3:end,:);
% P=[Q;P(3:end,:)];
M={};

for i=1:2:length(P)-2
    p=P([i i+1 i+2],:);
    q=Q([i i+1 i+2],:);
    MM=mobius_from_3points(p,q);
    M{end+1}=MobiusTrans(MM(1),MM(2),MM(3),MM(4));
    assert(norm(q-M{end}.map(p))<1e-10);
end
assert(length(FullP)==k*2-2);

temp=length(FullP)+2-(2:(k-1));
matching_vertices=[1 temp k (k-1):-1:2]';
% figure(1);
% clf
% hold on
% plot(FullP(:,1),FullP(:,2));
% plot(av(1),av(2),'bX');
% plot(s(1),s(2),'gX');
% plot(e(1),e(2),'rX');
%moving polygon so symmetric
dir=e-s;
if ~symmetric
    Mr=disc_mobius_translation(0,0);%(s(1)+1i*s(2),0);
else
Mr=disc_mobius_translation(s(1)+1i*s(2),0);
end
e2=Mr.map(e);
theta=-atan2(e(2),e(1))+pi/2;
if ~symmetric
    M2=disc_mobius_translation(0,-0);
M3=disc_mobius_translation(0/2,0);
else
M2=disc_mobius_translation(0,-theta);
e3=M2.map(e2);
M3=disc_mobius_translation(e3/2,0);
end


Mr=M3.compose(M2.compose(Mr));
FullP=Mr.map(FullP);
for i=1:length(M)
    M{i}=Mr.compose(M{i}.compose(Mr.inverse()));
end
% if nargout<4
%     return;
% end
[I,J]=meshgrid(1:length(FullP),1:length(FullP));
x=FullP(:,1);
y=FullP(:,2);
a=x(I);
b=y(I);
c=x(J);
d=y(J);
delta=2*((a-c).^2+(b-d).^2)./((1-a.^2-b.^2).*(1-c.^2-d.^2));
            D=acosh(1+delta);
            assert(all(all(abs(D-D')<1e-10)));
            dd=min(D,D(matching_vertices,:));
distance_matrix=dd(1:k,1:k);

end

