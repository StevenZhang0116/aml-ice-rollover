%--------------------------------------------------------------------------
% Edge tracking function
% Alternative Choice: OriginLab (https://www.originlab.com/)
% This function is quite messy and should work properly -- do not change unless have a
% strong reason... 
%
% Zihan Zhang, Courant Institute
% Updated June 2023
%--------------------------------------------------------------------------

function [grey,ashape,originalbw3,orig,ex,ey,outermostx,outermosty,polarx,polary,ltrap,rtrap,...
    lasttrap,prevarea,prevan,bestalpha_r,testx,testy,siga,edgethre] ...
    = func_track(v,f,timeInterval,cropped_interval,wst,wsb,left_adj,right_adj,threstrap,areaopen, ...
    waterareaopen, alphanum,lasttrap,prevarea,fr)

    frame = read(v,f); %import file?
    frame = imrotate(frame, 90); %rotate image to be vertical
    % crop the image to more focus on the flopping ice & around area
    croppedframe = imcrop(frame,cropped_interval); %what is cropped interval
    % find edges of ice
    % preprocessing
    grey = rgb2gray(croppedframe); %change to gray (why)
    grey = imadjust(grey); %?
    grey = imsharpen(grey); %?
    orig = grey;

    %% preliminary analysis
    % sensitivity threshold
    edgethre = 0.1; %how is this chosen / what does it do
    bw = edge(grey,'sobel',edgethre); %why 'sobel' (says for horizontal and vertical)
    sizeofImage = size(grey); 
    waterRange = [0,wst,sizeofImage(2),wsb-wst]; %how we choose parameters wsb, wst

    % delete pixels for some top/bottom side of image that is clearly not
    % part of the boundary detection
    % Takes cropped frame (grey) and removes points from bw (boundary
    % detection) that are 50 pixels(?) from the boundary
    bound1 = 50; %this bound chosen empirically? 
    for i=1:sizeofImage(2)
        for j=1:sizeofImage(1)
            if j<=bound1 || j>=sizeofImage(1)-bound1
                bw(j,i)=0;
            end
        end
    end

    waterpart = imcrop(grey, waterRange); %how waterRange chosen?

    %% Trap detection
    disp('Trap Detection/Iteration')
    % calculate the Gx, Gy gradient of cropped water surface area
    % [Gx,Gy] = imgradientxy(I) returns the directional gradients, Gx and Gy of the grayscale or binary image I.
    [Gx,~] = imgradientxy(waterpart); %in waterline, find grayscale gradients (large gradients ~ trap boundary probably)
    %specifically only compute Gx, which gives horizontal gradient at every
    %(kk,jj) in grid
    %ONLY horizontal gradients -- looking for horizontal cutoff of traps?
    %size(Gx) = size of waterpart (long x, thin y)
    gradthreshold = 100; %chosen empirically?
    cntupsize = []; cntdownsize = [];
    %Looking only within the water line region:
    for ii = 2:size(Gx,2)-1 
        %RESTART AND GO OVER THIN STRIP. Iterating over (x-waterpart)x3
        %array. Moving up over each ii...
        %Why do we iterate over thin strips of thickness 3???
        cntup = 0; cntdown = 0;
        for jj = ii-1:ii+1 %What are we counting? Don't we recount a lot of values this way...
            for kk=1:size(Gx,1) %iterate over all x values in waterline region
                if Gx(kk,jj) > gradthreshold
                    cntup = cntup + 1; %cntup = number of positive gradients > threshold
                end
                if Gx(kk,jj) < -gradthreshold
                    cntdown = cntdown + 1; %cntdown = number of |negative gradients| > threshold
                end
            end
        end
        %cntupsize(idx) = number of gradients over gradthreshold in
        %(x-waterpart)x3 array.
        cntupsize(end+1) = cntup; cntdownsize(end+1) = cntdown;
    end

    % adjustment of misdetection if the border line of trap is located in wrong part of the image
    % left trap
    [~,indup] = max(cntupsize); %% gives the index of the maximum cntupsize...
    % So identifies the horizontal strip of size (x-waterpart)x3 array with
    % the most number of big horizontal gradients.
    while indup >= 1/2*size(Gx,2) %Eliminate if indup is in the top half of the y-values. Why?
        cntupsize(indup) = 0;
        [~,indup] = max(cntupsize);
    end

    % possibly, dark large shadow of ice will be misidentified as trap; correction
    tt = 0; 
    while tt == 0
        %Isolate peaks in a certain region -- why this region
        %waterpart(20,1:indup)?
        %0.5?
        [~,ind] = findpeaks(im2single(waterpart(20,1:indup)),'MinPeakHeight',0.5);
        if isempty(ind) == 1 % if empty
            tt = 1; %if no peaks identified, skip
        else
            if mean(ind) < indup %Why does this mean bad left trap detection?
                disp("Bad LEFT Trap Detection")
                cntupsize(indup) = 0; %remove indup and try again.
                [~,indup] = max(cntupsize);
            end
        end
    end
    
    % right trap
    %Same thing but cntdown instead of cntup because white-->black instead
    %of black-->white.
    [~,inddown] = max(cntdownsize);
    while inddown <= 1/2*size(Gx,2)
        cntdownsize(inddown) = 0;
        [~,inddown] = max(cntdownsize);
    end
    
    %Why is this in here a second time?????
    %Shouldn't it be changed to be inddown, cntdownsize, right trap, etc?
    tt = 0; 
    while tt == 0
        [~,ind] = findpeaks(im2single(waterpart(20,indup:end)),'MinPeakHeight',0.5);
        if isempty(ind) == 0
            tt = 1;
        else
            if mean(ind) < indup
                disp("Bad LEFT Trap Detection")
                cntupsize(indup) = 0;
                [~,indup] = max(cntupsize);
            end
        end
    end

    % x coordinate of two vertical lines of trap
    % upcood: LHS of trap; downcood: RHS of trap
    origup = 1 + indup; origdown = 1 + inddown; 
    ltrap = origup + left_adj; %What are left_adj and right_adj?
    rtrap = origdown - right_adj;
    %ltrap and rtrap give indices for the y-range over which the trap is??

    % rollback system if misidentification occurs
    %What does this do?
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
    %Basically delete out gradient points that are identified as "trap"
    for i=1:size(Gx,2)
        for j=1:size(Gx,1)
            if i <= ltrap || i >= rtrap
                Gx(j,i)=0;
            end
        end
    end

    % find edge around water surface area
    %What does this do???
    waterpartedge = edge(Gx); 
    adjw = bwareaopen(waterpartedge,waterareaopen);
    %bwareaopen: Matlab function that creates new image...doing what?
    %waterareaopen?

    %Confused on this whole bit oops
    % connect edges??
    bw2 = bwmorph(bw, 'bridge');
    bw2 = bwmorph(bw2, 'diag');
    % delete unwanted edges??
    bw3 = bwareaopen(bw, areaopen);
    
    % stick together??
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
    [outermostx,outermosty] = fr.angle_sort(outermostx,outermosty);

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
    [testx,testy] = fr.angle_sort(testx,testy);


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


end




