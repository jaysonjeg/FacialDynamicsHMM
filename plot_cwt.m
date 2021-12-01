function [ax1,ax2] = plot_cwt(cmap_sig,values,period,Fs,coi,signif)
%{
Adapted from wt.m by Grinsted
Plot single CWT image array)

INPUTS:
cmap_sig is colormap, .e.g 'jet'
values is single CWT image array (nfreqs x ntimepoints)
period lists periods in seconds (1 x nfreqs)
Fs is sampling frequency (Hz)
coi (optional) has cone of influence
sig is optional array for significance contour plotting, where 1
indicates significant, 0 indicates nonsignificant. If sig is given, nonsig
values are greyscaled on an overlaid axis
%}

ax2=0; ax1=0;

t=(0:size(values,2)-1)/Fs;
plottype='Freq'; %'Freq' or 'Period'
logvalues=false; %Plot log of CWT values, or just plain values
alter_clim=false; %force change colour limits
alter_colorbarTicks=false;

if logvalues
    drawvalues=log2(values);
else
    drawvalues=values;
end

if nargin>=6
    ax1=axes;
else
    ax1=gca;
end

im1=imagesc(ax1,t,log2(period),drawvalues);   

if nargin>=6
    im1.AlphaData=signif;
end

Yticks = (2.^(fix(log2(min(period))):fix(log2(max(period)))));
Yticks=Yticks(1:2:end); %display less y values
if alter_clim
    clim=get(gca,'clim'); %center color limits around log2(1)=0
    clim=[-1 1]*max(clim(2),3);
    set(gca,'clim',clim)
end
%Set x/y labels
if strcmp(plottype,'Freq')
    YTickLabel=num2str(1./Yticks',1);
    %ylabel('Frequency (Hz)');
elseif strcmp(plottype,'Period')
    YTickLabel=num2str(Yticks',1);
    ylabel('Period (s)');
end
set(gca,'YLim',log2([min(period),max(period)]), ...
    'YDir','reverse', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',YTickLabel, ...
    'layer','top',...
    'xlim',[0,size(values,2)/Fs]);
%xlabel('Time (s)')

hold on

%cmap_sig=brewermap(100,'RdPu');
%cmap_sig=turbo;
%cmap_sig=inferno(100);


cmap_nonsig=rgb2gray(cmap_sig);
colormap(ax1,cmap_sig);


%Draw significant contours
if nargin>=6
    hold all;
    ax2=axes;    
    im2=imagesc(ax2,t,log2(period),drawvalues);
    im2.AlphaData=~signif;  
    linkaxes([ax1,ax2]);
    ax2.Visible='off';
    ax2.XTick=[]; ax2.YTick=[];
    colormap(ax2,cmap_nonsig);
    %cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]); 
    %cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]); 
     %{
    [c,h] = contour(t,log2(period),sig,[1 1],'k'); %#ok
    set(h,'linewidth',2)
    %}
end

%Draw COI
if nargin>=5
    hold on
    dt=1/Fs;
    %t2=(1:length(t))';
    t2=t';
    tt=[t2([1 1])-dt*.5;t2;t2([end end])+dt*.5];
    hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]),'w');
    set(hcoi,'alphadatamapping','direct','facealpha',.5)
    
end

%Set colorbar ticks
%HCB=colorbar;
if alter_colorbarTicks
    set(HCB,'ytick',-7:7);
    barylbls=rats(2.^(get(HCB,'ytick')'));
    barylbls([1 end],:)=' ';
    barylbls(:,all(barylbls==' ',1))=[];
    set(HCB,'yticklabel',barylbls);
end


hold off
set(gca,'box','on','layer','top');

end

