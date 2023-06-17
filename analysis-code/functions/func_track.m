%--------------------------------------------------------------------------
% Edge tracking function
% Alternative Choice: OriginLab (https://www.originlab.com/)
%
% Steven Zhang, Courant Institute
% Updated June 2023
%--------------------------------------------------------------------------

function [grey,ashape,originalbw3,orig,ex,ey,outermostx,outermosty,polarx,polary,ltrap,rtrap,...
    lasttrap,prevarea,prevan,bestalpha_r,testx,testy,siga,edgethre] ...
    = func_track(v,f,timeInterval,cropped_interval,wst,wsb,left_adj,right_adj,threstrap,areaopen, ...
    waterareaopen, alphanum,lasttrap,prevarea,fr)

    frame = read(v,f);
    frame = imrotate(frame, 90);
    % crop the image to more focus on the flopping ice & around area
    croppedframe = imcrop(frame,cropped_interval);
    % find edges of ice
    % preprocessing
    grey = rgb2gray(croppedframe);
    grey = imadjust(grey);
    grey = imsharpen(grey);
    orig = grey;

    %% preliminary analysis
    % sensitivity threshold
    edgethre = 0.1; 
    bw = edge(grey,'sobel',edgethre);
    sizeofImage = size(grey); 
    waterRange = [0,wst,sizeofImage(2),wsb-wst];

    % delete pixels for some top/bottom side of image that is clearly not
    % part of the boundary detection
    bound1 = 50;
    for i=1:sizeofImage(2)
        for j=1:sizeofImage(1)
            if j<=bound1 || j>=sizeofImage(1)-bound1
                bw(j,i)=0;
            end
        end
    end

    waterpart = imcrop(grey, waterRange);

    % calculate the Gx, Gy gradient of cropped water surface area
    [Gx,~] = imgradientxy(waterpart);
    gradthreshold = 100;
    cntupsize = []; cntdownsize = [];
    for ii = 2:size(Gx,2)-1
        cntup = 0; cntdown = 0;
        for jj = ii-1:ii+1
            for kk=1:size(Gx,1)
                if Gx(kk,jj) > gradthreshold
                    cntup = cntup + 1;
                end
                if Gx(kk,jj) < -gradthreshold
                    cntdown = cntdown + 1;
                end
            end
        end
        cntupsize(end+1) = cntup; cntdownsize(end+1) = cntdown;
    end

    % adjustment of misdetection if the border line of trap is located
    % in wrong part of the image
    [~,indup] = max(cntupsize);
    while indup >= 1/2*size(Gx,2)
        cntupsize(indup) = 0;
        [~,indup] = max(cntupsize);
    end
    [~,inddown] = max(cntdownsize);
    while inddown <= 1/2*size(Gx,2)
        cntdownsize(inddown) = 0;
        [~,inddown] = max(cntdownsize);
    end

    % x coordinate of two vertical lines of trap
    % upcood: LHS of trap; downcood: RHS of trap
    origup = 1 + indup; origdown = 1 + inddown; 
    ltrap = origup + left_adj;
    rtrap = origdown - right_adj;

    % rollback system if misidentification occurs
    if f == timeInterval(1)
        lasttrap = [ltrap,rtrap];
    else
        thistrap = [ltrap,rtrap];
        ratioset = thistrap./lasttrap;
        for i = 1:2
            if abs(ratioset(i)-1) > threstrap
                thistrap(i) = lasttrap(i);
            end
        end
        ltrap = thistrap(1);
        rtrap = thistrap(2);
        lasttrap = thistrap;
    end
    
    % delete the points in the box
    for i=1:size(Gx,2)
        for j=1:size(Gx,1)
            if i <= ltrap || i >= rtrap
                Gx(j,i)=0;
            end
        end
    end

    % find edge around water surface area
    waterpartedge = edge(Gx); 
    adjw = bwareaopen(waterpartedge,waterareaopen);

    % connect edges
    bw2 = bwmorph(bw, 'bridge');
    bw2 = bwmorph(bw2, 'diag');
    % delete unwanted edges
    bw3 = bwareaopen(bw, areaopen);
    
    % stick together
    for i=1:size(adjw,1)
        for j=1:size(adjw, 2)
            bw3(wst+i,j) = adjw(i,j);
        end
    end

    %% "raw" alphashape
    [prevarea,prevan,edges,bestalpha_r,ashape] = ...
        ashape_iteration(f,bw,bw3,prevarea,alphanum,areaopen,timeInterval,adjw,wst);

    ex = edges(:,1);
    ex(end+1) = ex(1);
    ey = edges(:,2);
    ey(end+1) = ey(1);
    
    % store original bw3
    originalbw3 = bw3;
    sizebw3 = size(bw3);

    %% outermost x 
    outermostbw3 = bw3;
    for k=1:sizebw3(1)
        therow = outermostbw3(k,:);
        nonzero = find(therow ~=0);
        % nonzero element
        sizenonzero = size(nonzero);
        % delete all other points except the leftmost and rightmost
        if sizenonzero(2)>2
            for ii=2:sizenonzero(2)-1
                outermostbw3(k,nonzero(ii)) = 0;
            end
        end
    end
    [outermosty,outermostx] = find(outermostbw3 ~= 0);

    % sort it using the angle wrt the centroid point
    xCenter = mean(outermostx); yCenter = mean(outermosty);
    angles = atan2d(outermosty-yCenter,outermostx-xCenter);
    [~,sai1] = sort(angles);
    outermostx = outermostx(sai1);
    outermosty = outermosty(sai1);

    %% test case: filter through x+y
    outermostbw3yy = bw3;
    for k=1:sizebw3(2)
        therow = outermostbw3yy(:,k);
        nonzero = find(therow ~=0);
        % nonzero element
        sizenonzero = size(nonzero);
        % delete all other points except the leftmost and rightmost
        if sizenonzero(1)>2
            for ii=2:sizenonzero(1)-1
                outermostbw3yy(nonzero(ii),k) = 0;
            end
        end
    end
    [testy,testx] = find(outermostbw3yy ~= 0);

    comm = [[testy;outermosty],[testx;outermostx]];
    % eliminate duplicated points 
    % comm = unique(comm,'row');

    shp = alphaShape(comm);
    siga = criticalAlpha(shp,'all-points');
    shp2 = alphaShape(comm,siga,'HoleThreshold',500000000);
    [~,edges] = boundaryFacets(shp2);
    testx = edges(:,2);
    testx(end+1) = testx(1);
    testy = edges(:,1);
    testy(end+1) = testy(1);
    
    % sort it using the angle wrt the centroid point
    xCenter = mean(testx); yCenter = mean(testy);
    angles = atan2d(testy-yCenter,testx-xCenter);
    [~,sai1] = sort(angles);
    testx = testx(sai1);
    testy = testy(sai1);


    %% polar 
    polarbw3 = originalbw3;
    [idy,idx] = find(polarbw3 ~= 0);
    cvh = convhull(idx,-idy);
    convx = idx(cvh); convy = -idy(cvh);
    [geo1,~,~] = fr.polygeom(convx,convy);
    % find the centroid use the convex hull
    xcm1 = geo1(2); ycm1 = geo1(3);
    % standardize the points coordinates
    outpx = idx - xcm1; 
    outpy = -idy - ycm1;
    % transfer to polar coordinates
    thetad = atan2d(outpy,outpx);
    rd = sqrt(outpx.^2+outpy.^2);
    % sort by theta
    [thetad,indx] = sort(thetad);
    rd = rd(indx);
    % rounding
    thetad = round(thetad);
    % pick the largest radius with respect to a certain angle
    angleset = []; radiusset = [];
    samplerate = 1;
    aang = -180:samplerate:180;
    for i = 1:length(aang)
        ang = aang(i);
        sett = (thetad == ang);
        radset = rd(sett);
        if size(radset,1) ~= 0
            maxrad = max(radset);
            angleset(end+1) = ang;
            radiusset(end+1) = maxrad;
        end
    end
    % change back to the x-y coordinate
    polarx = cosd(angleset).*radiusset + xcm1;
    polary = sind(angleset).*radiusset + ycm1;
    polary = -polary;

    % =================================== % 
    % interpolation and write r as a function of theta
    
    % Use Poly Interpolation (for now) 
%     rtind = 'P';
%     ind = 0; 
% 
%     % even distributed points for interpolation
%     aa = linspace(-180,180,natpara.numpt);
% 
%     if rtind(1) == 'F' && ind == 1
%         [nr,~,~,~] = intp.interpft(radiusset',angleset',10,natpara.numpt);
%         polarx_e = -cosd(aa).*nr' + xcm1;
%         polary_e = sind(aa).*nr' - ycm1;
%         figure()
%         plot(polarx_e,-polary_e)
%         hold on
%         plot(polarx,-polary)
%         sgtitle(f)
%         axis equal
%         close
%         pause;
% 
%     elseif rtind(1) == 'P' && ind == 1
%         % transform the units
%         % degree to radian
%         angleset = angleset./180.*pi; 
%         aa = aa./180.*pi;
%         % pixel to cm
%         radiusset = radiusset./natpara.rr; 
%         polarx = polarx./natpara.rr; polary = polary./natpara.rr;
%         xcm1 = xcm1/natpara.rr; ycm1 = ycm1/natpara.rr;
% 
%         % iterate through diff poly degree
%        [finaldegree,~] = intp.polyfit_inter(natpara.vrr,angleset,...
%            radiusset,aa,0);
%         
%         [ptheta,Sx] = polyfit(angleset,radiusset,finaldegree);
%         % px is the poly function of theta to describe radius
%         [finalr,~] = polyval(ptheta,aa,Sx);
%         
%         % transform back to xy coordinate
%         polarx_e = cos(aa).*finalr + xcm1;
%         polary_e = sin(aa).*finalr + ycm1;
%         polary = -polary; 
% 
%         % calculate curvature using r-theta interpolation
%         % http://mathonline.wikidot.com/the-curvature-of-plane-polar-curves
%         f_1de = polyder(ptheta);
%         f_2de = polyder(f_1de);
% 
%         param1 = fr.polyadd(2*(conv(f_1de,f_1de)),conv(ptheta,ptheta));
%         param1 = fr.polyadd(param1,-conv(ptheta,f_2de));
%         param2 = fr.polyadd(conv(f_1de,f_1de),conv(ptheta,ptheta));
%         kappa = abs(polyval(param1,aa)) ./ ((abs(polyval(param2,aa))).^(3/2));
% 
%         [peaks,loc] = findpeaks(kappa,'MinPeakDistance',20);
%         cre = (peaks > mean(kappa) + 1*std(kappa));
%         peaks = peaks(cre); loc = loc(cre);
% 
%         figure('units','normalized','outerposition',[0 0 1 1]);
% 
%         % plot
%         subplot(1,3,1)
%         plot(angleset,radiusset)
%         hold on
%         plot(aa,finalr)
%         title(['degree=',num2str(finaldegree),'    R-\theta interpolation'])
%         legend('real','intp');
%         xlabel('angle (radian)')
%         ylabel('radius (cm)')
%         hold off
% 
%         subplot(1,3,2)
%         plot(aa,kappa);
%         hold on
%         for i = 1:length(loc)
%             theloc = loc(i);
%         end
%         xlabel('angle (radian)')
%         title('\kappa in R-\theta calculation')
% 
%         subplot(1,3,3)
%         hold on
%         scatter(polarx_e,polary_e,'red')
%         scatter(polarx,polary,'green')
%         scatter(polarx_e(1),polary_e(1),'filled','blue')
%         scatter(polarx_e(10),polary_e(10),'filled','yellow')
%         scatter(xcm1,ycm1,'*')
%         for i = 1:length(loc)
%             theloc = loc(i);
%         end
%         title(['Polar Interpolation with ',num2str(length(polarx_e)),...
%             ' points'])
%         legend('interpolation','polar outermost','start pt','ind pt',...
%             'center')
%         axis equal
%         xlabel('x axis (cm)')
%         ylabel('y axis (cm)')
% 
%         % mean radius
%         mr = mean(finalr);
%         [xxx,yyy] = fr.mkcircle(xcm1,ycm1,mr,bkt);
%         scatter(xxx,yyy);
% 
%         hold off
%         
%         sgtitle(f)
%         close;
%     end
end




