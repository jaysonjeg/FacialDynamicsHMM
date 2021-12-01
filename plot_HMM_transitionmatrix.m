function pplot=plot_HMM_transitionmatrix(adj,layout,prc,thiscolormap)
adj=network_cutoffpercentile(adj,prc);
G=digraph(adj,'omitselfloops'); %Px1 for default, Px2 for top 20 percentile of connections
LWidths = 8*G.Edges.Weight/max(G.Edges.Weight);
if strcmp(layout,'layered')
    pplot=plot(G,'Layout','layered',...
    'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths,...
    'ArrowSize',15,'ArrowPosition',0.4,'MarkerSize',8,'NodeFontSize',20);
elseif strcmp(layout,'force')
    pplot=plot(G,'Layout','force','WeightEffect','inverse',...
        'LineWidth',LWidths,...
        'ArrowSize',15,'ArrowPosition',0.5,'MarkerSize',11,'NodeFontSize',15,...
        'NodeColor',thiscolormap);
    
    %'EdgeLabel',G.Edges.Weight);
end
pplot.EdgeLabel=cellfun(@(x) sprintf('%.2f',eval(x)),pplot.EdgeLabel,'UniformOutput',false);

%Move node label positions
text(pplot.XData+0.13, pplot.YData+0.08 ,pplot.NodeLabel, ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment', 'left',...
    'FontSize', 15)
pplot.NodeLabel = {}; 

end

