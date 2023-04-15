function []=mystimulusannotate(colours,thickness,only_comedy)
%Called by Figure 3 in melancholia_analysis_mat2cube

if nargin==2
    only_comedy=false;
end

ratio=0.5; %0 to 1. 1 means under the axis. 0 means over.
y0=min(get(gca,'ylim'));
[xaf1,yaf1] = ds2nfu(0,y0);
[xaf2,yaf2]=ds2nfu(125.6,y0);
[xaf3,yaf3]=ds2nfu(125.6+168.84,y0);
[xaf4,yaf4]=ds2nfu(352.1200,y0);
annotation('rectangle',[xaf1,yaf1-thickness*ratio,xaf2-xaf1,thickness],'Color','none','FaceColor',colours(1,:));

if ~only_comedy
    annotation('rectangle',[xaf2,yaf1-thickness*ratio,xaf3-xaf2,thickness],'Color','none','FaceColor',colours(2,:));
    annotation('rectangle',[xaf3,yaf1-thickness*ratio,xaf4-xaf3,thickness],'Color','none','FaceColor',colours(3,:));
end

end

