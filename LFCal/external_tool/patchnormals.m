function N = patchnormals(FV) 
%Vertex normals of a triangulated mesh, area weighted, left-hand-rule 
% N = patchnormals(FV) -struct with fields, faces Nx3 and vertices Mx3 
%N: vertex normals as Mx3
%face corners index 
% added by Vincent

A = FV.Faces(:,1); 
B = FV.Faces(:,2); 
C = FV.Faces(:,3);
%face normals 
n = cross(FV.Vertices(A,:)-FV.Vertices(B,:),FV.Vertices(C,:)-FV.Vertices(A,:)); %area weighted
%vertice normals 
N = zeros(size(FV.Vertices)); %init vertix normals 
for i = 1:size(FV.Faces,1) %step through faces (a vertex can be reference any number of times) 
N(A(i),:) = N(A(i),:)+n(i,:); %sum face normals 
N(B(i),:) = N(B(i),:)+n(i,:); 
N(C(i),:) = N(C(i),:)+n(i,:); 
end