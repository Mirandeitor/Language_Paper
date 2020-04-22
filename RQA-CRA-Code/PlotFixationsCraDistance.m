function PlotFixationsCraDistance(fixation1,fixation2,radius,imageFilename)

    f=figure(1);
    clf

    % Read image
    img=imread(imageFilename);
    img=128+floor(img/2);
    imagesc(img);
    hold on

    nfix=size(fixation1,1);

    % Plot fixations
    plot(fixation1(:,1),fixation1(:,2),'k');
    plot(fixation2(:,1),fixation2(:,2),'r');
    axis ij; %puts y-coor at top left (to fix for how eyelink does it)
    axis equal;
    axis([0 1024 0 768]); %set to image size
    set(gca,'XTick',[]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTick',[]);
    set(gca,'YTickLabel',[]);
    hold on

    % Plot circles
    ang=0:0.05:2*pi; 
    xp=radius*cos(ang);
    yp=radius*sin(ang);

    for i = 1:nfix-1
        for j=i+1:nfix
        dist = sqrt((fixation1(i,1)-fixation2(j,1))^2+(fixation1(i,2)-fixation2(j,2))^2);
        if dist < radius
            plot(fixation1(i,1)+xp, fixation1(i,2)+yp,'k')
            hold on
        end
    end

end
