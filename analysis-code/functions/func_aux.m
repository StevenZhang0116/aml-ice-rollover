%--------------------------------------------------------------------------
% All registered auxiliary functions for miscellaneous purposes in ice-rollover experiments
%
% Steven Zhang, Courant Institute
% Updated July 2023
%--------------------------------------------------------------------------

% functions registeration
function fr = func_aux
    fr.sigmoid = @sigmoid;
    fr.polygeom = @polygeom; 
    fr.sort_by_angle = @sort_by_angle;
    fr.peak_select = @peak_select;
    fr.polyadd = @polyadd;
    fr.detect_watersurface = @detect_watersurface;
    fr.calibration_watersurface = @calibration_watersurface;
    fr.calibration_video = @calibration_video;
    fr.mkcircle = @mkcircle;
    fr.incircle = @incircle;
    fr.minboundcircle = @minboundcircle;
    fr.p_poly_dist = @p_poly_dist;
    fr.circlefit = @circlefit;
    fr.mkngon = @mkngon;
    fr.matrix_rotation = @matrix_rotation;
end

% create data points around the regular pentagon
function [xx,yy,origx,origy,LL] = mkngon(x,y,r,bkt,intp,nn,pert,pertr,pertm,kappa,inwardflag)
    % only one type of perturbation could occur
    assert(kappa * pert == 0);

    ex = []; ey = []; 
    vang = []; % set of angles of vertex
    for i = 1:nn+1
        ex(end+1) = x + r*cos(i*2*pi/nn);
        ey(end+1) = y + r*sin(i*2*pi/nn);
        vang(end+1) = i*2*pi/nn;
    end

    sideset = sqrt(diff(ex).^2+diff(ey).^2);
    sidelength = sideset(1);

    ss = sqrt((diff(ex).^2)+(diff(ey).^2));
    LL = ss(1) * nn;
    pt = intp.interparc(bkt+1,ex,ey,'linear');
    dinc = (360-0)/bkt; % increment of angle degree

    % harmonic perturbations
    if pert ~= 0 
        xx = pt(:,1);
        yy = pt(:,2);
        r = sqrt(xx.^2 + yy.^2); % length=bkt
        thetad = atan2d(yy,xx);
        % sort by degree
        [thetad,indx] = sort(thetad);
        r = r(indx); % thetad --- r

        % must be continuous (harmonic formulation)
        ind1 = (thetad >= pertr(1));
        part1 = thetad(ind1);
        ind2 = (part1 <= pertr(2));
        part2 = part1(ind2);
        
        nor_pertr = part2/(part2(end)-part2(1)) * 2 * pi; 
        part_rinc = sin(nor_pertr) * pertm; % part2 --- part_rinc

        for i = 1:length(thetad)
            ang = thetad(i);
            ff = find(part2 == ang);
            if ff ~= 0
                r(i) = r(i) + part_rinc(ff);
            end
        end

        xx = r.*cosd(thetad);
        yy = r.*sind(thetad);
    else
        xx = pt(:,1);
        yy = pt(:,2);
    end

    origx = xx;
    origy = yy;

    % use fix curvature perturbation
    if kappa ~= 0
        makerr = 1/kappa;
        spt = 16384;
        gridstep = 2*pi/spt;

        newxx = [];
        newyy = [];
    
        for jjj = 1:length(ex)-1
            [xxcc,yycc] = mkcircle(x,y,makerr,spt); % circle with given curvature
            startpt = [xxcc(1),yycc(1)];
            reldistset = [];
            for hhh = 1:length(xxcc)
                distt = sqrt((xxcc(hhh)-startpt(1)).^2+(yycc(hhh)-startpt(2)).^2);
                reldistset(end+1) = distt;
            end
            reldistset = abs(reldistset-sidelength);
            % chrod length equal to the polygon side
            tol = 1e-4;
            indd = reldistset < tol;
            idx = find(indd,1,'first');
            % generate (bkt+1)/nn points during this arc
            frac = round((bkt+1)/nn);
            allpt = spt * (frac)/idx;
            [xxcc,yycc] = mkcircle(x,y,makerr,allpt);
            tracsec = [xxcc(1:frac);yycc(1:frac)]';

            % check for each side
            disp(['side:',num2str(jjj)])
            move = [ex(jjj)-tracsec(end,1),ey(jjj)-tracsec(end,2)];
            tracsec = tracsec + move;
            midpp = length(tracsec)/2;

            rotangle = 0:0.001:2*pi;
            if inwardflag == 0 % outward
                if jjj == 1
                    nextjjj = nn;
                else
                    nextjjj = jjj-1;
                end
            else % inward
                if jjj == nn
                    nextjjj = 1;
                else
                    nextjjj = jjj+1;
                end
            end
            distset = [];
            for kkk = 1:length(rotangle)
                rota = rotangle(kkk);
                [rotx,roty] = matrix_rotation(ex(jjj),ey(jjj),tracsec(:,1),tracsec(:,2),rota);
                dist = sqrt((rotx(1)-ex(nextjjj)).^2 + (roty(1)-ey(nextjjj)).^2);
                distset(end+1) = dist;
            end
            [~,opt] = min(distset);
            optangle = rotangle(opt);
            [rotx,roty] = matrix_rotation(ex(jjj),ey(jjj),tracsec(:,1),tracsec(:,2),optangle);

            newxx = [newxx;rotx];
            newyy = [newyy;roty];
        end

        % illustration
        ill = 1;
        if ill == 1
            figure()
            axis equal
            hold on
            scatter(newxx,newyy);
            plot(xx,yy,'color','r','LineStyle','--')
            legend('Fixed curvature','Original Polygon')
        end

        % replace the points
        xx = newxx;
        yy = newyy;
        
    end
end

% calibration analysis2
% focus on the water surface analysis
% return the coordinates of leftmost & rightmost points in the top side
% of water surface (i.e. the intersection with the trap), and the water 
% thickness
function [head2top] = calibration_watersurface(videoName,interval,f)
%     v = VideoReader(videoName);
%     frame = read(v,f);
% 
%     grey = imrotate(frame,90);
%     grey = imcrop(grey,interval);
%     
%     % conduct the manual calibration; 
%     figure()
%     set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
%     title('calibration')
%     imshow(grey)
%     hold on
%     line([0,1080],[50,50],'LineWidth',2)
% 
%     h = imdistline;
%     fcn = makeConstrainToRectFcn('imline',get(gca,'XLim'),get(gca,'YLim'));
%     setDragConstraintFcn(h,fcn);   
%     
%     disp('ALL UNITS IN PIXELS');

%     head2top = input('distance between top water surface to the top of the cropped frame:');
    head2top = 183.49;
end

function [rotx,roty] = matrix_rotation(x,y,xset,yset,ang)
    stdx = xset-x;
    stdy = yset-y;
    rotx = stdx*cos(ang)-stdy*sin(ang);
    roty = stdx*sin(ang)+stdy*cos(ang);
    rotx = rotx+x;
    roty = roty+y;
end

% create data points around the circle 
function [xx,yy] = mkcircle(x,y,r,bkt)
    ang = 0:2*pi/bkt:2*pi; 
    xp = r*cos(ang);
    yp = r*sin(ang);
    xx = x+xp;
    yy = y+yp;
end

% smoothed step function
function val = sigmoid(z,zL,zR,H,ds)
    val = (zL - zR)./(1 + exp((z - H)/ds)) + zR;
end

% least square circle fitting
function [xo,yo,R] = circlefit(x,y)
    x = x(:);
    y = y(:);
    A = [-2*x -2*y ones(length(x),1)];
    x = A\-(x.^2+y.^2);
    xo=x(1);
    yo=x(2);
    R = sqrt(xo.^2+yo.^2-x(3));
end

% calculate properties of polygon
function [geom, iner, cpmo] = polygeom(x, y) 
    %   POLYGEOM( X, Y ) returns area, X centroid,
    %   Y centroid and perimeter for the planar polygon
    %   specified by vertices in vectors X and Y.
    %
    %   [ GEOM, INER, CPMO ] = POLYGEOM( X, Y ) returns
    %   area, centroid, perimeter and area moments of 
    %   inertia for the polygon.
    %   GEOM = [ area   X_cen  Y_cen  perimeter ]
    %   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
    %     u,v are centroidal axes parallel to x,y axes.
    %   CPMO = [ I1     ang1   I2     ang2   J ]
    %     I1,I2 are centroidal principal moments about axes
    %         at angles ang1,ang2.
    %     ang1 and ang2 are in radians.
    %     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv
     
    if ~isequal(size(x), size(y))
        error('X and Y must be the same size');
    end
     
    % temporarily shift data to mean of vertices for improved accuracy
    xm = mean(x);
    ym = mean(y);
    x = x - xm;
    y = y - ym;
      
    % summations for CCW boundary
    xp = x([2:end 1]);
    yp = y([2:end 1]);
    a = x.*yp - xp.*y;
     
    A = sum(a) /2;
    xc = sum((x+xp).*a) /6/A;
    yc = sum((y+yp).*a) /6/A;
    Ix = sum((y.*y+y.*yp+yp.*yp).*a)/12;
    Iy = sum((x.*x+x.*xp+xp.*xp).*a)/12;
    Ixy = sum((x.*yp+2*x.*y+2*xp.*yp+xp.*y).*a) /24;
     
    dx = xp - x;
    dy = yp - y;
    P = sum(sqrt(dx.*dx+dy.*dy));
     
    % check for CCW versus CW boundary
    if A < 0
        A = -A;
        Ix = -Ix;
        Iy = -Iy;
        Ixy = -Ixy;
    end
     
    % centroidal moments
    Iu = Ix-A*yc*yc;
    Iv = Iy-A*xc*xc;
    Iuv = Ixy-A*xc*yc;
    J = Iu+Iv;
    
    % replace mean of vertices
    x_cen = xc+xm;
    y_cen = yc+ym;
    Ix = Iu+A*y_cen*y_cen;
    Iy = Iv+A*x_cen*x_cen;
    Ixy = Iuv+A*x_cen*y_cen;
     
    % principal moments and orientation
    I = [ Iu  -Iuv ;
         -Iuv   Iv ];
    [ eig_vec, eig_val ] = eig(I);
    I1 = eig_val(1,1);
    I2 = eig_val(2,2);
    ang1 = atan2( eig_vec(2,1), eig_vec(1,1) );
    ang2 = atan2( eig_vec(2,2), eig_vec(1,2) );
     
    % return values
    geom = [ A  x_cen  y_cen  P ];
    iner = [ Ix  Iy  Ixy  Iu  Iv  Iuv ];
    cpmo = [ I1  ang1  I2  ang2  J ];
end



% sort set of points by rotation angle (thres)
function [newpt,indret] = sort_by_angle(pt,thres)
    xx = pt(:,1);
    yy = pt(:,2);
    rr = sqrt(xx.^2 + yy.^2);
    aa = atan2d(yy,xx);
    for i=1:length(aa)
        if aa(i) < 0
            aa(i) = aa(i) + 360;
        end
    end
    [~,sa] = sort(aa);
    rr = rr(sa); aa = aa(sa);
    [~,lc] = min(aa < thres);
    
    xxx = cosd(aa).*rr;
    yyy = sind(aa).*rr;

    fpx = xxx(1:lc-1); fpy = yyy(1:lc-1);
    xxx = [xxx(lc:end)]; yyy = [yyy(lc:end)];
    xxx = [xxx;fpx]; yyy = [yyy;fpy];

    newpt(:,1) = xxx;
    newpt(:,2) = yyy;

    % sorting index
    stpt = [xxx(1);yyy(1)];
    distostpt = sqrt((xx-stpt(1)).^2+(yy-stpt(2)).^2);
    [~,ctt] = min(distostpt);
    indret = [ctt:1:length(xx),1:1:ctt-1];
   
end

% select peaks of given data based on certain criteria
function [peaks,loc] = peak_select(aset,pd,pp)
    [peaks,loc] = findpeaks(aset,'MinPeakDistance',pd,...
            'MinPeakProminence',pp);
    % what if the start/end point around is the peak (global maximum) 
    % where the peak function is uncapable to identify
    if aset(1) == aset(end)
        kk = [];
        if aset(end) > aset(end-1)
            kk(1,:) = aset(1:20);
        end
        if(aset(1) > aset(2))
            kk(2,:) = aset(end-19:end);
        end
        ind = [0,0];
        for i=1:size(kk,1)
            [~,indd] = max(kk(i,:));
            ind(i) = indd;
        end
        if(ind(1) == 1 && ind(2) == length(aset))
            ind(2) = [];
        end
        ind = ind(ind ~= 0);
        for i=1:size(ind,2)
            peaks(end+1) = aset(ind(i));
            loc(end+1) = ind(i);
        end
    end

    peaks = unique(peaks);
    loc = unique(loc);
end

% calibration analysis; transform pixel unit to centimeter
function [ratio] = calibration_video(videoName,f)
    v = VideoReader(videoName);
    frame = read(v,f);
    
    % conduct the manual calibration; 
    figure()
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    title('calibration')
    imshow(frame)
    h = imdistline;
    fcn = makeConstrainToRectFcn('imline',get(gca,'XLim'),get(gca,'YLim'));
    setDragConstraintFcn(h,fcn);   
    
    disp('ALL UNITS IN PIXELS');
    distinpixel = input('value from read: ');
    distincm = input('centimeter used: ');
    
    ratio = distinpixel/distincm;
    ratio = ratio / 10^(-2); % change unit to m
    disp(['The conversion ratio is: ', num2str(ratio)]);
end

% addition of polynomial functions
function [result] = polyadd(p1,p2)
    result = [zeros(1, size(p1,2)-size(p2,2)) p2] + [zeros(1, size(p2,2)-size(p1,2)) p1];
end

% detailed analysis of water surface and trap range and set this part of 
% region to be 0 in the matrix 
% [outdated]
function [val1,val2,direction_ind] = detect_watersurface( ...
    videoreader,waterSurfaceTop,waterSurfaceBot,crp,numofFrame)
    frame = read(videoreader,numofFrame,"native");
    frame = imrotate(frame, 90);
    croppedframe = imcrop(frame,crp);

    % preprocessing of image
    grey = rgb2gray(croppedframe);

    waterrange = [0 waterSurfaceTop 1080 waterSurfaceBot-waterSurfaceTop];
    waterpart = imcrop(grey,waterrange);

    waterpart = imadjust(waterpart);
    waterpart = imlocalbrighten(waterpart,'AlphaBlend',true);
    waterpart = imsharpen(waterpart);
    waterpart = imreducehaze(waterpart);

    % Contrast-limited adaptive histogram equalization (CLAHE)
    waterpart = adapthisteq(waterpart,'clipLimit',0.02,'Distribution',...
        'rayleigh');
    % Enhance contrast using histogram equalization
    waterpart = histeq(waterpart);

    % calculate gradient
    [Gmag, Gdir] = imgradient(waterpart,'sobel');

    bd1 = edge(waterpart,'prewitt','nothinning');
    bd1 = bwmorph(bd1, 'bridge');
    bd1 = bwmorph(bd1, 'diag');
    bd1 = bwareaopen(bd1, 20);

    figure(1);    
    imshow(waterpart)

    % hough transform to detect the horizontal line
    % since the detection is very unstable, we choose to use the longest 
    % line as the primal reference
    [H,theta,rho] = hough(bd1);
    P = houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))));
    lines = houghlines(bd1,theta,rho,P,'FillGap',20,'MinLength',50);
    hold on
    max_len = 0;
    yheight = [];
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       y1 = lines(k).point1(2);
       y2 = lines(k).point2(2);

       % only care about "parallel" lines
       if abs(y1-y2) <= 5
           yheight(end+1)=y1;
           yheight(end+1)=y2;
        
           % Determine the endpoints of the longest line segment
           len = norm(lines(k).point1 - lines(k).point2);
           if (len > max_len)
              max_len = len;
              xy_long = xy;
           end
       end
    end
    plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');

    val1 = xy_long(:,1)+crp(1);
    val2 = xy_long(:,2)+crp(2);
    direction_ind = 0; % set the direction index

    % detect whether the baseline is located around the upper or lower
    % water surface. Draw the rectangle therefrom. 
    % -1 meaning upward expansion, 1 meaning downward expansion
    if mean(val2) > (waterSurfaceBot-waterSurfaceTop)/2
        direction_ind = -1;
    else
        direction_ind = 1;
    end

    hold off
    close;
end

% calculation incircle
% compute the maximal in-circle of the polygonal convex hull of a set of 
% points in the plane
function [C,R] = incircle(x,y)
    if (nargin < 1) || (nargin > 2)
      error('INCIRCLE:improperarguments', ...
        'incircle requires exactly 1 or 2 arguments')
    elseif (nargin == 2) && isvector(x) && isvector(y) && (numel(x) ...
            == numel(y))
      % a pair of vectors
      xy = [x(:),y(:)];
      % compute the hull. I prefer convhulln.
      edges = convhulln(xy);
    elseif (nargin == 1) && (size(x,2) == 2)
      % a single list of points as rows of x
      xy = x;
      edges = convhulln(xy);
    elseif (nargin == 2) && (size(x,2) == 2) && (size(y,2) == 2)
      % y must be a list of edges from convhulln
      xy = x;
      edges = y;
    elseif (nargin == 2) && (size(x,2) == 2) && isvector(y) && (y(1) ...
            == y(end))
      % y must be a list of edges from convhull
      xy = x;
      I = y(:);
      % in case the first point was not wrapped in the polygon
      if I(1) ~= I(end)
        I = [I;I(1)];
      end
      edges = [I(1:(end-1)),I(2:end)];
    else
      % none of the forms I allow seem to fit.
      % give up and throw an error.
      error('INCIRCLE:invaliddata', ...
        ['x and y do not seem to fit into any of the allowed forms for ' ...
        'incircle input'])
    end
    ne = size(edges,1);
    
    % the end points of each edge are...
    A = xy(edges(:,1),:);
    B = xy(edges(:,2),:);
    
    % the normal vector to each edge
    N = (B - A)*[0 1;-1 0];
    
    % normalize to unit length
    L = sqrt(sum(N.^2,2));
    N = N./[L,L];
    
    % a central point inside the hull itself
    C0 = mean(A,1);
    
    % ensure the normals point inwards. While
    % I can use the signed area for a hull as
    % generated by convhull (suggestion by Bruno)
    % this test is also efficient, and it deals
    % with the case of a hull provided from some
    % unknown and less trusted source.
    k = sum(N.*bsxfun(@minus,C0,A),2) < 0;
    N(k,:) = -N(k,:);
    
    % equality constraints defined by the dot products
    % with the normal vectors.
    Aeq = [N,zeros(ne,1),-eye(ne)];
    beq = sum(N.*A,2);
    
    % lower bounds only for the slack variables
    LB = [-inf; -inf; 0; zeros(ne,1)];
    % there are no upper bounds
    UB = [];
    
    % inequalities defined by the slack variable
    % constraints
    A = [zeros(ne,2),ones(ne,1),-eye(ne)];
    b = zeros(ne,1);
    
    % the objective is just -R
    f = zeros(ne+3,1);
    f(3) = -1;
    
    % just turn off the output message
    options = optimset('linprog');
    options.Display = 'off';
    
    % linprog does all of the hard work.
    result = linprog(f,A,b,Aeq,beq,LB,UB,[],options);
    
    % unpack the circle parameters
    C = result(1:2)';
    R = result(3);

end

% calculation circumscribed circle
% compute the minimum radius enclosing circle of a set of (x,y) pairs
function [center,radius] = minboundcircle(x,y,hullflag)
    % default for hullflag
    if (nargin<3) || isempty(hullflag)
      hullflag = true;
    elseif ~islogical(hullflag) && ~ismember(hullflag,[0 1])
      error 'hullflag must be true or false if provided'
    end
    
    % preprocess data
    x=x(:);
    y=y(:);
    
    % not many error checks to worry about
    n = length(x);
    if n~=length(y)
      error 'x and y must be the same sizes'
    end
    
    % start out with the convex hull of the points to
    % reduce the problem dramatically. Note that any
    % points in the interior of the convex hull are
    % never needed.
    if hullflag && (n>3)
      edges = convhulln([x,y]);
    
      % list of the unique points on the convex hull itself
      % convhulln returns them as edges
      edges = unique(edges(:));
    
      % exclude those points inside the hull as not relevant
      x = x(edges);
      y = y(edges);
        
    end
    
    % now we must find the enclosing circle of those that
    % remain.
    n = length(x);
    
    % special case small numbers of points. If we trip any
    % of these cases, then we are done, so return.
    switch n
      case 0
        % empty begets empty
        center = [];
        radius = [];
        return
      case 1
        % with one point, the center has radius zero
        center = [x,y];
        radius = 0;
        return
      case 2
        % only two points. center is at the midpoint
        center = [mean(x),mean(y)];
        radius = norm([x(1),y(1)] - center);
        return
      case 3
        % exactly 3 points
        [center,radius] = enc3(x,y);
        return
    end
    
    % more than 3 points.
    
    % Use an active set strategy.
    aset = 1:3; % arbitrary, but quite adequate
    iset = 4:n;
    
    % pick a tolerance
    tol = 10*eps*(max(abs(mean(x) - x)) + max(abs(mean(y) - y)));
    
    % Keep a list of old sets as tried to trap any cycles. we don't need to
    % retain a huge list of sets, but only a few of the best ones. Any cycle
    % must hit one of these sets. Really, I could have used a smaller list,
    % but this is a small enough size that who cares? Almost always we will
    % never even fill up this list anyway.
    old.sets = NaN(10,3);
    old.rads = inf(10,1);
    old.centers = NaN(10,2);
    
    flag = true;
    while flag
      % have we seen this set before? If so, then we have entered a cycle
      aset = sort(aset);
      if ismember(aset,old.sets,'rows')
        % we have seen it before, so trap out
        center = old.centers(1,:);
        radius = old.radius(1);
        
        % just reset flag then continue, and the while loop will terminate
        flag = false;
        continue
      end
      
      % get the enclosing circle for the current set
      [center,radius] = enc3(x(aset),y(aset));
      
      % is this better than something from the retained sets?
      if radius < old.rads(end)
        old.sets(end,:) = sort(aset);
        old.rads(end) = radius;
        old.centers(end,:) = center;
            
        % sort them in increasing order of the circle radii
        [old.rads,tags] = sort(old.rads,'ascend');
        old.sets = old.sets(tags,:);
        old.centers = old.centers(tags,:);
      end
      
      % are all the inactive set points inside the circle?
      r = sqrt((x(iset) - center(1)).^2 + (y(iset) - center(2)).^2);
      [rmax,k] = max(r);
      if (rmax - radius) <= tol
        % the active set enclosing circle also enclosed
        % all of the inactive points.
        flag = false;
      else
        % it must be true that we can replace one member of aset
        % with iset(k). Which one?
        s1 = [aset([2 3]),iset(k)];
        [c1,r1] = enc3(x(s1),y(s1));
        if (norm(c1 - [x(aset(1)),y(aset(1))]) <= r1)
          center = c1;
          radius = r1;
          
          % update the active/inactive sets
          swap = aset(1);
          aset = [iset(k),aset([2 3])];
          iset(k) = swap;
          
          % bounce out to the while loop
          continue
        end
        s1 = [aset([1 3]),iset(k)];
        [c1,r1] = enc3(x(s1),y(s1));
        if (norm(c1 - [x(aset(2)),y(aset(2))]) <= r1)
          center = c1;
          radius = r1;
          
          % update the active/inactive sets
          swap = aset(2);
          aset = [iset(k),aset([1 3])];
          iset(k) = swap;
          
          % bounce out to the while loop
          continue
        end
        s1 = [aset([1 2]),iset(k)];
        [c1,r1] = enc3(x(s1),y(s1));
        if (norm(c1 - [x(aset(3)),y(aset(3))]) <= r1)
          center = c1;
          radius = r1;
          
          % update the active/inactive sets
          swap = aset(3);
          aset = [iset(k),aset([1 2])];
          iset(k) = swap;
          
          % bounce out to the while loop
          continue
        end
        
        % if we get through to this point, then something went wrong.
        % Active set problem. Increase tol, then try again.
        tol = 2*tol;
        
      end
      
    end
    
    % =======================================
    %  begin subfunction
    % =======================================
    function [center,radius] = enc3(X,Y)
    % minimum radius enclosing circle for exactly 3 points
    %
    % x, y are 3x1 vectors
    
    % convert to complex
    xy = X + sqrt(-1)*Y;
    
    % just in case the points are collinear or nearly so, get
    % the interpoint distances, and test the farthest pair
    % to see if they work.
    Dij = @(XY,i,j) abs(XY(i) - XY(j));
    D12 = Dij(xy,1,2);
    D13 = Dij(xy,1,3);
    D23 = Dij(xy,2,3);
    
    % Find the most distant pair. Test if their circumcircle
    % also encloses the third point.
    if (D12>=D13) && (D12>=D23)
      center = (xy(1) + xy(2))/2;
      radius = D12/2;
      if abs(center - xy(3)) <= radius
        center = [real(center),imag(center)];
        return
      end
    elseif (D13>=D12) && (D13>=D23)
      center = (xy(1) + xy(3))/2;
      radius = D13/2;
      if abs(center - xy(2)) <= radius
        center = [real(center),imag(center)];
        return
      end
    elseif (D23>=D12) && (D23>=D13)
      center = (xy(2) + xy(3))/2;
      radius = D23/2;
      if abs(center - xy(1)) <= radius
        center = [real(center),imag(center)];
        return
      end
    end
    
    % if we drop down to here, then the points cannot
    % be collinear, so the resulting 2x2 linear system
    % of equations will not be singular.
    A = 2*[X(2)-X(1), Y(2)-Y(1); X(3)-X(1), Y(3)-Y(1)];
    rhs = [X(2)^2 - X(1)^2 + Y(2)^2 - Y(1)^2; ...
           X(3)^2 - X(1)^2 + Y(3)^2 - Y(1)^2];
         
    center = (A\rhs)';
    radius = norm(center - [X(1),Y(1)]);
    
    end

end

function [d,x_poly,y_poly] = p_poly_dist(x, y, xv, yv) 
    % If (xv,yv) is not closed, close it.
    xv = xv(:);
    yv = yv(:);
    Nv = length(xv);
    if ((xv(1) ~= xv(Nv)) || (yv(1) ~= yv(Nv)))
        xv = [xv ; xv(1)];
        yv = [yv ; yv(1)];
    %     Nv = Nv + 1;
    end
    % linear parameters of segments that connect the vertices
    % Ax + By + C = 0
    A = -diff(yv);
    B =  diff(xv);
    C = yv(2:end).*xv(1:end-1) - xv(2:end).*yv(1:end-1);
    % find the projection of point (x,y) on each rib
    AB = 1./(A.^2 + B.^2);
    vv = (A*x+B*y+C);
    xp = x - (A.*AB).*vv;
    yp = y - (B.*AB).*vv;
    % Test for the case where a polygon rib is 
    % either horizontal or vertical. From Eric Schmitz
    id = find(diff(xv)==0);
    xp(id)=xv(id);
    clear id
    id = find(diff(yv)==0);
    yp(id)=yv(id);
    % find all cases where projected point is inside the segment
    idx_x = (((xp>=xv(1:end-1)) & (xp<=xv(2:end))) | ((xp>=xv(2:end)) & ...
        (xp<=xv(1:end-1))));
    idx_y = (((yp>=yv(1:end-1)) & (yp<=yv(2:end))) | ((yp>=yv(2:end)) & ...
        (yp<=yv(1:end-1))));
    idx = idx_x & idx_y;
    % distance from point (x,y) to the vertices
    dv = sqrt((xv(1:end-1)-x).^2 + (yv(1:end-1)-y).^2);
    if(~any(idx)) % all projections are outside of polygon ribs
       [d,I] = min(dv);
       x_poly = xv(I);
       y_poly = yv(I);
    else
       % distance from point (x,y) to the projection on ribs
       dp = sqrt((xp(idx)-x).^2 + (yp(idx)-y).^2);
       [min_dv,I1] = min(dv);
       [min_dp,I2] = min(dp);
       [d,I] = min([min_dv min_dp]);
       if I==1 %the closest point is one of the vertices
           x_poly = xv(I1);
           y_poly = yv(I1);
       elseif I==2 %the closest point is one of the projections
           idxs = find(idx);
           x_poly = xp(idxs(I2));
           y_poly = yp(idxs(I2));
       end
    end
    if(inpolygon(x, y, xv, yv)) 
       d = -d;
    end
end






