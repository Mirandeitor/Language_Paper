function PlotFixationsDistance(xfix,yfix,radius,imageFilename,plotFilename)

    f=figure(1);
    clf

    % Read image
    img=imread(imageFilename);
    img=128+floor(img/2);
    imagesc(img);
    hold on

    nfix=length(xfix);

    % Plot fixations
    pl = plot(xfix,yfix,'k');
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
        %plot a circle when the distance between two successive fixations i and
        %j are less than the specified radius (i.e., when they recur)
        for j = i+1:nfix
            dist = sqrt(( xfix(i)-xfix(j) )^2 + ( yfix(i)-yfix(j) )^2);
            %if dist < radius
                plot(xfix(i)+xp, yfix(i)+yp,'k')
                hold on
            %end
        end
    end

    fullFilename=[plotFilename '.png'];
    print(f,'-dpng',fullFilename);

end
