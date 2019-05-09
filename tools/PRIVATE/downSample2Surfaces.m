function cortexN = downSample2Surfaces(Nodes1, Faces1, Nodes2, Faces2, r1, r2)
cortex1.vertices = Nodes1;
cortex1.faces    = Faces1;
cortexN1         = reducepatch(cortex1,r1);

cortex2.vertices = Nodes2;
cortex2.faces    = Faces2;
cortexN2         = reducepatch(cortex2, r2);
cortexN.vertices = [cortexN1.vertices;cortexN2.vertices];
cortexN.faces = [cortexN1.faces;cortexN2.faces+size(cortexN1.vertices,1)];
end