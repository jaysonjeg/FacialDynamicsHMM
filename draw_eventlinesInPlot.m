function [] = draw_eventlinesInPlot(event_times,framerate)
hold on
for i=1:length(event_times)
    xline(event_times(i)*framerate,'LineWidth',2);
end
hold off
end