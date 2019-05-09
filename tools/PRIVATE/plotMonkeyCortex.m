function [h]= plotMonkeyCortex(cortex, inds, j, conf)
conf.ifFig=1;conf.TITLE='Sources near the electrodes';
patchVis(cortex, cortex.Vertices(inds(j),:), ones(1,length(j)), conf)
hold on
text(cortex.Vertices(inds(j),1), cortex.Vertices(inds(j),2), cortex.Vertices(inds(j),3), num2str([1:1:125]'),'FontSize',10,'Color', 'k');
view(-90,0); 
cameratoolbar;
h =gca;
end