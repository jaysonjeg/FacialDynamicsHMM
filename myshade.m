function []=myshade(t,ylow,yhigh,colour,FaceAlpha)
%In a plot with x-values given by 't', shade values between 'ylow' and
%'yhigh'

if all(size(ylow)==[1,1]) %if single numbers are given
    ylow=repmat(ylow,1,length(t));
    yhigh=repmat(yhigh,1,length(t));
end

%remove NaNs at beginning and end of ylow and yhigh
first_index=1; last_index=length(ylow);
if isnan(ylow(1))
    first_index=find(~isnan(ylow),1);
end
if isnan(ylow(end))
    last_index=length(ylow)+1-find(~isnan(flip(ylow)),1);
end
ylow=ylow(first_index:last_index);
yhigh=yhigh(first_index:last_index);
t=t(first_index:last_index);

if nargin<5
    FaceAlpha=0.3;
end
x2 = [t, fliplr(t)];
inBetween = [yhigh, fliplr(ylow)];
hold on;
fill(x2, inBetween, colour,'FaceAlpha',FaceAlpha,'EdgeColor','none'); %or 'patch'
hold off
end