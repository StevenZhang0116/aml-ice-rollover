%--------------------------------------------------------------------------
% Body Dynamics calculation
% Mainly focusing on center of mass, inertia, buoyancy force, and torque
% Calculations are well-aligned to Scott Weady's simulation
% Outdated due to no clear tendency is observed; replaced by pure simulation
%
% Zihan Zhang, Courant Institute
% Updated June 2023
%--------------------------------------------------------------------------

function bdyn = func_bodydynamics 
    bdyn.calbd_integral = @calbd_integral;
    bdyn.calbody_integral = @calbody_integral;
end

% use boundary integral to calculate the physical quantities
function [crosscnt,equi,totaltad,modtorset,moveset,Hret,xcm,ycm] = calbd_integral(pt,finalL,orig,...
    ex,ey,cx,cy,natpara,fr,ind)
    disp('==BOUNDARY==')
    global rho_water rho_air rho_ice g
    
    startdeg = 0;
    enddeg = 360-ind;
    totaltad = startdeg:ind:enddeg;
    moveset = []; % vertical movement (to resolve instability)
    modtorset = []; % modified torque set
    % fourier modes
    iks = 2*pi*1i * [(0:natpara.numpt/2),(-natpara.numpt/2+1:-1)]'; 
    x = pt(:,1); y = pt(:,2); % interpolation points
    
    H = -natpara.wcs; % water surface level, by calibration
    Ns = length(x);
    s = linspace(0,1,Ns+1)'; 
    s = s(1:end-1); 
    ds = s(2) - s(1); %arclength

    figure()
    subplot(1,2,1)
    imshow(orig)
    hold on
    scatter(ex * natpara.rr,ey * natpara.rr,5,'r','filled'); axis equal
    hold off
    title('Back Image and Alpha Boundary')
    subplot(1,2,2)
    plot(x,y);axis equal;
    hold on; scatter(cx,cy); hold off
    title('Interpolation Result')
    close
    
    ltt = length(totaltad);
    % store the coordinate of iceberg after certain rotation
    rrcell = cell(1,ltt);

    for i = 1:ltt
        thetad = totaltad(i);

        xs = real(ifft(iks.*fft(x)));
        ys = real(ifft(iks.*fft(y)));
        % == calculate the area == %
        A2 = abs(sum(x.*ys)*ds);
        % == calculate the center of mass, only calculate once when theta=0 == %
        if thetad == 0
            xcm = 0.5/A2*sum((x.^2).*ys)*ds; 
            ycm = -0.5/A2*sum((y.^2).*xs)*ds;
        end

        % rotation
        opx = x-xcm; opy = y-ycm; 
        % multiply with the rotation matrix
        px = opx*cosd(thetad)-opy*sind(thetad);
        py = opx*sind(thetad)+opy*cosd(thetad);
        % == coordinates after rotation == %
        rx = px+xcm; ry = py+ycm; 

        rrcell{i} = [rx,ry];

        % == calculate buoyancy force, where the center of mass does not matter == %
        rxx = rx; ryy = ry;
        ys = @(dz) (real(ifft(iks.*fft(ryy+dz))));
        if thetad == 0
            smoother2 = @(dz) (fr.sigmoid(ryy,rho_water,rho_air,H+dz,ds));
        else
            smoother2 = @(dz) (fr.sigmoid(ryy+dz,rho_water,rho_air,H,ds));
        end
        
        Fb = @(dz) (-g * sum(rxx.*ys(dz).*(rho_ice - smoother2(dz)))*ds);
        % find smallest buoyancy force (≈0)
        absoFb = @(dz) (abs(Fb(dz)));
        
        optset = optimset('TolX',1e-10);
        [param,~] = fminsearch(absoFb,0,optset);

        if thetad == 0
            % modify the initial water line height (theta=0) st the buoyance force=0
            H = H + param; 
            Hret = round(H,5);
            param = 0;               
        end

        moveset(end+1) = param; % optimal lift distance

        % == calculate torque == %
        rxx = rx-xcm; ryy = ry + param; % shift vertically
        ys = real(ifft(iks.*fft(ryy)));
        smoother1 = fr.sigmoid(ryy,rho_water,rho_air,H,ds);
        T = -g/2*ds*sum((rho_ice - smoother1).*(rxx.^2).*ys);

        % == calculate modified torque == %
        modtorque = T;
        modtorset(end+1) = modtorque;
    end

    % modify the graph presentation
    mdpt = ltt/2;
    totaltad = totaltad - 180;
    moveset = [moveset(mdpt+1:end),moveset(1:mdpt)];
    modtorset = [modtorset(mdpt+1:end),modtorset(1:mdpt)];
    
    % Find equilibrium points
    equi = [];
    for i = 1:length(modtorset)-1
        if modtorset(i) * modtorset(i+1) < 0
            equi(end+1) = i;
        end
    end

    % how many roots of "mod" torque
    crosscnt = length(equi);

end

% calculate the boundary dynamic using body integral of all grid points
% inside the alphaShape
function [crosscnt,equi,totaltad,modtorset,moveset,Hret,xcm,ycm,oldxcm] = calbody_integral(orig,...
    ashape,natpara,pt,flg,ind,appendheight,inter)
    %

    %
    disp('==BODY==')
    global rho_water rho_air rho_ice g 
    [~,aedges] = boundaryFacets(ashape);
    ax = aedges(:,1); ay = aedges(:,2); % alpha boundary
    if flg ~= 0
        ax = pt(:,1); ay = pt(:,2);
    end
    % the grid size is very fine, so that all vertex are located on the
    % boundary
    xl = min(ax):inter:max(ax);
    yl = min(ay):inter:max(ay);
    gridpoints = [];
    for i = 1:length(xl)
        for j = 1:length(yl)
            gridpoints(end+1,1) = xl(i); gridpoints(end,2) = yl(j);
        end
    end

    % find grid points inside the alphashape
    if flg == 0
        tf = inShape(ashape,gridpoints(:,1),gridpoints(:,2));
        ispt = gridpoints(tf,:);
    elseif flg ~= 0
        in = inpolygon(gridpoints(:,1),gridpoints(:,2),ax,ay);
        ispt = gridpoints(in,:);
    end

    dix = max(ispt(:,1))-min(ispt(:,1));
    diy = max(ispt(:,2))-min(ispt(:,2));

    ispt = ispt./natpara.rr;

    figure()
    scatter(ispt(:,1),ispt(:,2))
    axis equal
    close

    % dx: uniform grid in 2 directions
    % =========
    % =========
    dx = ((1/natpara.rr) * (1/(1/inter)))^2; 

    H = -natpara.wcs; % water surface level, by calibration

    startdeg = 0;
    enddeg = 360-ind;
    totaltad = startdeg:ind:enddeg;
    ltt = length(totaltad);

    moveset = [];
    modtorset = [];
    pointloccell = cell(1,ltt);
    
    for i = 1:ltt
        thetad = totaltad(i);

        % calculate center of mass
        % all rotation is around the center of mass
        if thetad == 0
            xcm = mean(ispt(:,1));
            oldxcm = xcm;
            ycm = mean(ispt(:,2));
        end

        % rotation
        opx = ispt(:,1)-xcm; opy = ispt(:,2)-ycm; 
        % multiply with the rotation matrix
        px = opx*cosd(thetad)-opy*sind(thetad);
        py = opx*sind(thetad)+opy*cosd(thetad);
        % == coordinates after rotation == %
        rx = px+xcm; ry = py+ycm; 

        pointloccell{i} = [rx,ry];

        % calculate buoyancy
        % density function
        if thetad == 0
            dens = @(dz) (ry <= H+dz) * rho_water + (ry > H+dz) * rho_air; 
        else 
            dens = @(dz) (ry+dz <= H) * rho_water + (ry+dz > H) * rho_air; 
        end

        % buoyancy force calculation, Equation 4.4
        Fb = @(dz) sum((rho_ice - dens(dz)) * -g * dx);

        absoFb = @(dz) (abs(Fb(dz)));
        optset = optimset('TolX',1e-10);
        [param,~] = fminsearch(absoFb,0,optset);
        
        if thetad == 0
            H = H + param;
            % use the water height calculated from boundary integral (more accurate)? 
            H = appendheight; 
            Hret = round(H,5);
            param = 0;
        end

        moveset(end+1) = param;

        % optimal shift
        optdz = param;

        % update density function using the updated waterheight
        % not “pointer”
        if thetad == 0
            dens = @(dz) (ry <= H+dz) * rho_water + (ry > H+dz) * rho_air; 
        end
        
        % find "more stable" center of mass
        if thetad == 0
            dista = @(modxcm) (rx-modxcm);
            torq = @(modxcm) sum((dista(modxcm) * -g).*(rho_ice-dens(0)) * dx);
            abstorq = @(modxcm) abs(torq(modxcm));
            optset = optimset('TolX',1e-10);
            [movex,~] = fminsearch(abstorq,xcm,optset);
            % update center of mass
            xcm = movex;
            modtorset(end+1) = 0;
        else
            % calculate torque
            dista = rx-xcm;
            % Equation 4.5
            torq = sum((dista * -g).*(rho_ice-dens(optdz)) * dx);
            modtorset(end+1) = torq;
        end
    end

    % plot
    mdpt = ltt/2;
    totaltad = totaltad - 180;
    moveset = [moveset(mdpt+1:end),moveset(1:mdpt)];
    modtorset = [modtorset(mdpt+1:end),modtorset(1:mdpt)];

    equi = [];
    for i = 1:length(modtorset)-1
        if modtorset(i) * modtorset(i+1) <= 0
            if (isempty(equi) == 1) || (equi(end) ~= i-1)
                equi(end+1) = i;
            end
        end
    end

    crosscnt = length(equi);

    % number of equilibrium point
    showequi = 1;
    if showequi == 1
        nnn = length(equi);
        figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])
        for iiiii = 1:nnn
            aa = equi(iiiii);
            pp = totaltad(aa); % points after rotation
            hh = moveset(aa); % move distance
        
            opx = ispt(:,1)-xcm; opy = ispt(:,2)-ycm; 
            px = opx*cosd(pp)-opy*sind(pp);
            py = opx*sind(pp)+opy*cosd(pp);
            rx = px+xcm; ry = py+ycm;
        
            subplot(1,nnn,iiiii)
            scatter(rx,ry+hh);
            hold on
            yline(Hret);
            title(['\theta=',num2str(pp)],'FontSize',12)
            axis equal
            hold off
        end
    end
end
























