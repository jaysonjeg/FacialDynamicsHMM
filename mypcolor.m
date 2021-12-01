function out=mypcolor(x,y,myimage)
%{
pcolor plot with axis labels.
Add column/row to right and bottom so it doesn't cut off
Argument 'y' is optional
%}
if size(y,2)==1
    y=y'; %make into row vector
end
if nargin<3
    myimage=y;
    y=1:height(myimage);
end
myimage2=[myimage,myimage(:,end)];
myimage3=[myimage2;myimage2(1,:)];
out=pcolor([x,x(end)+1],[y,y(end)+1],myimage3);
shading flat; colormap('jet'); colorbar;
end