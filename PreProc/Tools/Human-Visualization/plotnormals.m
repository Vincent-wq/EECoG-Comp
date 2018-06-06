function [n,T] = plotnormals(fig,S,opt,color,size)

if fig
    figure;
    axis equal
end
switch opt
    case 'faces'
        trep = TriRep(S.faces,S.vertices);
        n = faceNormals(trep)*sign(size);
        T = (S.vertices(S.faces(:,1),:)+S.vertices(S.faces(:,2),:)+S.vertices(S.faces(:,3),:))/3;
    case 'vertices'
        n = getNormals(S)*sign(size);
        T = S.vertices;
    otherwise
        error('option must be ''faces'' or ''vertices''');
end
hold on; quiver3(T(:,1),T(:,2),T(:,3),n(:,1),n(:,2),n(:,3),abs(size),color);