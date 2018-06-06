function plotgrads(iffig,lighton,alpha,grad,R,r,colors,types)

G1 = grad(:,1:3);
G2 = grad(:,4:6);
if size(grad,2)==9
    N = grad(:,7:9);
else
    N = G2-G1;
    N = N./repmat(sqrt(dot(N,N,2)),1,3);
end
Nc = 5;
T = torus(R,r,Nc);
T1 = T;
T2 = T;
n = size(grad,1);
if iffig
    figure;
else
    hold on;
end
for i = 1:n
    D = D_Matrix(N(i,:));
    M1 = [D G1(i,:)'; 0 0 0 1];
    T1.vertices = applyM(T.vertices,M1);
    M2 = [D G2(i,:)'; 0 0 0 1];
    T2.vertices = applyM(T.vertices,M2);
    switch types
        case 'in'
            plotsurf(T1,0,alpha,0,'none',colors(1,:));
        case 'out'
            plotsurf(T2,0,alpha,0,'none',colors(2,:));
        case 'both'
            plotsurf(T1,0,alpha,0,'none',colors(1,:));
            plotsurf(T2,0,alpha,0,'none',colors(2,:));
    end
end
if lighton
    camlight(0,180)%,'infinite')
    camlight(0,0)%,'infinite')
    material dull
end

function D = D_Matrix(mu)

% Rotates the axis z to the vector mu 

u = mu/norm(mu);
ux = u(1);
uy = u(2);
uz = u(3);

if all(mu(:)==[0 0 1]')
    D = eye(3);
else
    m = sqrt(1-uz^2);
    D = [ux*uz/m -uy/m ux;...
        uy*uz/m  ux/m uy;...
        -m     0   uz];
end
    
function T = torus(R,r,Nc,Nc2)

dcr = 2*pi/Nc;
if nargin == 3
    dcR = r/R*dcr;
else
    dcR = 2*pi/Nc2;
end
[Cr,CR] = ndgrid(0:dcr:2*pi,0:dcR:2*pi);
T.faces =  delaunay(Cr,CR);
X = (R+r.*cos(Cr)).*cos(CR);
Y = (R+r.*cos(Cr)).*sin(CR);
Z = r.*sin(Cr);
T.vertices = [X(:) Y(:) Z(:)];
T = unique_surf(T,1e-3);