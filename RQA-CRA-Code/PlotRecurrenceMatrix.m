function PlotRecurrenceMatrix(recurrenceMatrix,plotFilename)
nSeq=max(size(recurrenceMatrix));

[x,y,~]=find(recurrenceMatrix);

h=figure(2);
clf;

plot(x,y,'o','MarkerSize',5.0,'MarkerEdgeColor','r','MarkerFaceColor','r');
xlim([1 nSeq]);
ylim([1 nSeq]);
axis square;
title(['Recurrence Plot '])
%title(['Recurrence Plot ' plotFilename])

%print(h,'-dpng',plotFilename);

end
