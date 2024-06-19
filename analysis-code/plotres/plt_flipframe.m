%--------------------------------------------------------------------------
% Plot frames just after a pervious flip and just before the next flip
%
% Zihan Zhang, Courant Institute
% Updated May 2023
%--------------------------------------------------------------------------

function [] = plt_flipframe(alltInv,pttracker,indInv,closeoverind,wl)
    if closeoverind == 1 
        close all 
    end
        
    time = alltInv;
    frame = pttracker;
    indicator = indInv; 
    
    assert(length(alltInv) == length(frame))

    nn = length(alltInv);
    [ha, pos] = tight_subplot(round(nn/2),2,[.01 .03],[.1 .01],[.01 .01]);
    for ii = 1:nn
        thistime = time{ii};
        thisind = indicator{ii};
        combind = [];
        for jj = 1:length(thisind)
            combind(end+1) = thisind{jj}(1); 
        end
        findvalid = find(combind==1);

        assert(length(findvalid) == 2)

        axes(ha(ii)); 
        pts = frame{ii};
        initpt = pts{findvalid(1)};
        endpt = pts{findvalid(2)};
        hold on
        plot(initpt(:,1),initpt(:,2),'LineStyle','--','LineWidth',2)
        plot(endpt(:,1),endpt(:,2),'LineWidth',2)
        yline(wl,'LineWidth',1.5,'LineStyle','-.')
        legend(num2str(thistime(findvalid(1))),num2str(thistime(findvalid(2))))
        axis equal
        hold off
    end
    set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
end

