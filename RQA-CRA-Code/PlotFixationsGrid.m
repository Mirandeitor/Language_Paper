function PlotFixationsGrid(xfix,yfix,gridsize,imageFilename,plotFilename)

f=figure(1);
clf

% Read image
img=imread(imageFilename);
img=128+floor(img/2);
imagesc(img);
hold on

%plot simple x,y fixations
pl = plot(xfix,yfix,'k');
axis ij; %puts y-coor at top left (to fix for how eyelink does it)
axis equal;
axis([0 1024 0 768]); %set to image size
set(gca,'XTick',0:gridsize:1024);
set(gca,'XTickLabel',[]);
set(gca,'YTick',0:gridsize:768);
set(gca,'YTickLabel',[]);
grid on
hold on

fullFilename=[plotFilename '.png'];
print(f,'-dpng',fullFilename);
