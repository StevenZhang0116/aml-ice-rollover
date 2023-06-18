%--------------------------------------------------------------------------
% Main function for ice flipping experiments 
% Boundary extraction and calculation of physical quantities
%
% Steven Zhang, Courant Institute
% Updated June 2023
%--------------------------------------------------------------------------


close all
clear
clc

%% Basic Setting
setting 
addpath(genpath('datas'));
addpath('plotres');
addpath('functions');
addpath('externs');
mkdir(['datas/flip-moments/',foldername])
disp('Create Folder for Record Flip-Around Moments')


%% Data Preloading
% load fliptime 
fliptimeload = 1;
if fliptimeload == 1
    date_experiment = foldername(1:end-1);
    flippath = 'data-rolltime/';
    adjoint = strcat(flippath,date_experiment,'.mat');
    flipdata = load(adjoint).result;
    % ignore flips after 23 minutes [roughly]
    flipdata = flipdata(flipdata<23*60);
    % also rollback several seconds
%     flipdata = flipdata-5; 
    % transform seconds to frame number
    flipdata = flipdata.*rfr;
    % add sf into it, as the indicator of initial shape
    fliptimestamp = [sf,flipdata];
    tInv = fliptimestamp;
end

% timestamp for calculating melting rate
meltrateind = 1;
if meltrateind == 1
    allmc = cell(1,length(tInv)-2);
    allnormalmc = cell(1,length(tInv)-2);
    pttracker = cell(1,length(tInv)-2);
    alltInv = cell(1,length(tInv)-2);
    indInv = cell(1,length(tInv)-2);
    centerlst = cell(1,length(tInv)-2);
    for i = 2:length(tInv)-1
        % the manually selected second is slighly before each flip
        cf = i;
        delay = 5;
        bf = tInv(cf)+2*delay*rfr; bf2 = tInv(cf)+3*delay*rfr; 
        af = tInv(cf+1)-delay*rfr; af2 = tInv(cf+1)-2*delay*rfr;
        inttflip = 5*rfr;
        % select "moment" between 2 flips
        flipInv = [bf:inttflip:bf2,af2:inttflip:af];
        alltInv{1,i-1} = flipInv;
        % indicator cell [maybe have miscellaneous usages]
        theind = cell(1,length(flipInv)); indInv{1,i-1} = theind; 
        % save overall tracking points
        ptt = cell(1,length(flipInv)); pttracker{1,i-1} = ptt;
        % measuring melting rate set that record the boundary underneath water
        mc = cell(1,length(flipInv)); allmc{1,i-1} = mc;
        % save normal vector
        normalmc = cell(1,length(flipInv)); allnormalmc{1,i-1} = normalmc;
        % save center of mass
        centerpt = cell(1,length(flipInv)); centerlst{1,i-1} = centerpt;
    end
end

%% start the main for loop to frames in the video
% change bounds in generating data with manual operations
for jj = 8:length(alltInv)
tInv = alltInv{jj};
tflipdiff = (tInv(end)-tInv(1))/rfr;
% specifiy the initial frame to the end
for iiii = 1:length(tInv) 
    f = tInv(iiii);
    disp(f/v.FrameRate)
    
    %% Edge tracking
    [grey,ashape,originalbw3,orig,ex,ey,outermostx,outermosty,polarx,polary,...
        upcood,downcood,lasttrap,prevarea,prevan,bestalpha_r,testx,testy,siga,edgethre] ...
        = func_track(v,f,tInv,cropped_interval,waterSurfaceTop,waterSurfaceBot, ...
        left_adj,right_adj,threstrap,areaopen,waterareaopen,alphanum,lasttrap,...
        prevarea,fr);

    % plot
    plt_grpind = 0;
    if plt_grpind == 1
        figure('units','normalized','outerposition',[0 0 1 1]);

        % normal
        subplot(1,4,1);
        imshow(orig);
        hold on;
        plot(ex,-ey,'r','LineWidth',2);
        title('Alpha Boundary Edge');
        alphaxy = [ex,ey];
        alphatrackcell{1,iiii} = alphaxy;
    
        subplot(1,4,2);
        imshow(originalbw3);
        title('Edges Detected')
        % store it
        edgetrackcell{1,iiii} = originalbw3;

        % OUTERMOST XY FORMULATION % 
        subplot(1,4,3)
        imshow(orig);
        hold on
        plot(outermostx,outermosty,'red','LineWidth',2)
        hold off
        title('XY outermost boundary')
        outermostxy = [outermostx,-outermosty];
        xytrackcell{1,iiii} = outermostxy;

        % OUTERMOST R \theta FORMULATION %
        subplot(1,4,4)
        imshow(orig);
        hold on
        plot(polarx*natpara.rr,-polary*natpara.rr,'red','LineWidth',2)
        hold off
        title('R\theta outermost boundary')
        polarxy = [polarx;-polary];
        rthetatrackcell{1,iiii} = polarxy;

        aa = f/(60*25);
        minutes = floor(aa);
        seconds = (aa-minutes)*60;
        sgtitle(["Frame" string(f),'Time:',num2str(minutes), ...
            num2str(seconds)]);

        close
    end

    if savefig == 1
        figure()
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        subplot(2,2,1)
        imshow(originalbw3)
        hold on
        line([upcood,upcood],[waterSurfaceTop,waterSurfaceBot],'LineWidth',2)
        line([downcood,downcood],[waterSurfaceTop,waterSurfaceBot],'LineWidth',2)
        title(edgethre)
        % raw alpha shape
        subplot(2,2,2)
        imshow(orig);
        hold on;
        scatter(ex,-ey,'r','LineWidth',2);
        line([upcood,upcood],[waterSurfaceTop,waterSurfaceBot],'LineWidth',2)
        line([downcood,downcood],[waterSurfaceTop,waterSurfaceBot],'LineWidth',2)
        title(["dummy alphashape: ", num2str(bestalpha_r)])  
        % raw outermost x
        subplot(2,2,3)
        imshow(orig);
        hold on
        scatter(outermostx,outermosty,'r','LineWidth',2)
        % arbitrary test case
        subplot(2,2,4)
        imshow(orig);
        hold on
        scatter(testx,testy,'r','LineWidth',2)
        title(["adaptive alphashape + X&Y filter", num2str(siga)])  
        % saveas(gcf,['./bd-tracking/',num2str(iiii),'.jpeg'])
        wantthis = input('Use this IMAGE or not? ');
        assert((wantthis==0) || (wantthis==1))
        close
    end

    %% manual editing of extracted boundary points
    if wantthis == 1 && manuelw == 1
        [newx,newy] = func_manual(orig,ex,-ey,fr);
    end

    %% Interpolation    
    % choose which type of detection we aim to use
    % ALPHASHAPE: ex, ey
    % XY: outermostx,outermosty
    % POLAR: polarx,polary
    if dec(1) == 'X'
        ex = outermostx; ey = outermosty; 
    elseif dec(1) == 'P'
        ex = polarx; ey = polary; 
        ex = ex'; ey = ey';
    % otherwise is the ALPHASHAPE option
    elseif dec(1) == 'A'
        ey = -ey;
    end

    expixel = ex; eypixel = ey;

    % use manuelly operated points to calculate
    if wantthis == 1 && manuelw == 1
        ex = newx; ey = newy; 
    end

    % calibration video rescale
    ex = ex/natpara.rr;
    ey = ey/natpara.rr;

    % Calculate the arclength set
    N = length(ex);
    eex = [ex;ex(1)];
    eey = [ey;ey(1)];

    ds = sqrt((eex(2:end)-eex(1:end-1)).^2+(eey(2:end)-eey(1:end-1)).^2);
    s = zeros(N,1); % arclength accumulation
    for n = 1:N-1
        s(n+1) = s(n) + ds(n);
    end
    finalL = s(end); % total arclength

    % equidistributed arclength
    xx = linspace(0,finalL,natpara.numpt+1); 
    % discard the final point, identical with the first one
    xx = xx(1:end-1); 

    % conduct interpolation based on the chosen method
    [pt,res] = func_allinterpolation(intp,natpara,ex,ey,s,xx);

    % unpack Fourier parameters
    if method(1) == 'F'
        kappa = res{1}; kappa_h = res{2}; curv = res{3}; x_e = res{4}; y_e = res{5}; 
        fd = res{6}; sd = res{7};
    end

    %% Calculate center of mass
    % pt is the interpolation points with x,y coordinates
    polyin = polyshape(pt);
    [cx,cy] = centroid(polyin); % center of mass -- 1

    [geom, iner, cpmo] = fr.polygeom(pt(:,1),pt(:,2));
    A1 = geom(1); % area of polygon
    areaSet(iiii) = A1;
    cx2 = geom(2); cy2 = geom(3); % center of mass-- 2
    inertia = iner(4)+iner(5); % calculate intertia -- 1

    % store the normalized interpolation point (-center of mass)
    intpcell{1,iiii} = [pt(:,1)-cx,pt(:,2)-cy]; 

        % save manully operated points to dataframe with timestamp
    if wantthis == 1 && manuelw == 1
        savethis = input('Save This File? ');
        if savethis == 1
            global rfr
            ts = f/rfr; % (s)
            of = ['datas/flip-moments/',foldername,num2str(ts),'.mat']; % output folder
            save(of,"pt","ts","cx","cy","ex","ey");
            disp(['Saved Data to ', of]); 
        end
    end

    %% Calculate harmonic moment (outdated)
    % x_h,y_h: fourier coefficients
    % x_e,y_e: fourier interpolation 

    krange = 0:1:20;
    harmind = 0;
    if harmind == 1
        [aa,mset] = func_harmonic(krange,cx,cy,x_e,y_e,xx,fd);
    end

    %% Calculate melting rate at specific frame
    meind = 1;
    if meind == 1
        % height of bottom side of water surface (cm)
        wtbt = (-waterSurfaceBot)/natpara.rr;
        
        % unpack
        ux = pt(:,1); uy = pt(:,2); 
        wl = -0.02995; % water line
        wlineset(end+1) = wl;

        % save tracking points
        pttracker{jj}{1,iiii} = pt-[cx,0];

        xs_e = fd(:,1); ys_e = fd(:,2);
        % points below the surface
        indd = (uy < wtbt);
        % normalize the position so that centered at (0,0)
        ux = ux(indd) - cx; uy = uy(indd) - cy; 
        xs_e = xs_e(indd); ys_e = ys_e(indd);
        % calculate unit normal vector
        nx = -ys_e ./ sqrt(xs_e.^2+ys_e.^2);
        ny = xs_e ./ sqrt(xs_e.^2+ys_e.^2);
        meltpt = [ux,uy];

        [newmeltpt,indret] = fr.sort_by_angle(meltpt,90);
        allmc{jj}{1,iiii} = newmeltpt; % store the frame
        indInv{jj}{1,iiii} = wantthis; 
        allnormalmc{jj}{1,iiii} = [nx(indret),ny(indret)]; % store the normal vector
        centerlst{jj}{1,iiii} = [cx,cy];
    end

    %% Calculate shape dynamics

    % for simplicity of calculation, traclonsform back the unit! 
    backpt = pt * natpara.rr / 100; 
    backcx = cx * natpara.rr / 100;
    backcy = cy * natpara.rr / 100;
    backex = ex * natpara.rr / 100;
    backey = ey * natpara.rr / 100;

    dynind = 0;
    if dynind == 1
        % =================================== %  
        % Approach 1
        % criteria of curvature check [only using the curvature value]
        % 1. Convex vertex
        % 2. Larger than mean+1*std
        % =================================== %  

        if backpt(1,:) ~= backpt(end,:)
            backpt = [backpt;backpt(1,:)];
        end
        [L2,R2,K2] = cf.curvature(backpt); % 3-pts method
        ccva = 1./R2;

        % directly using fourier coefficients
        % roughly same effect
        if method(1) == 'F'
            ccva = [curv;curv(1)]/mean(curv); 
        end
    
        [peaks,loc] = fr.peak_select(ccva,200,1.75);
    
        % statistics for curvature calculation
        stddiv = 1;
        meanc = mean(ccva);
        stdc = std(ccva);
    
        curve_threshold = 0;
    
        indlst = zeros(length(peaks)); % direction check
        stdlst = zeros(length(peaks)); % compartive value check 
        abslst = zeros(length(peaks)); % absolute value check 
        for i = 1:length(peaks)
            peakpt = [backpt(loc(i),1),backpt(loc(i),2)];
            vect = [K2(loc(i),1),K2(loc(i),2)];
            eps = 0.1;
            inpt = [peakpt(1)+vect(1)*eps, peakpt(2)+vect(2)*eps];
            in = inpolygon(inpt(:,1),inpt(:,2),backpt(:,1),backpt(:,2));
            indlst(i) = in;
            stdlst(i) = (ccva(loc(i)) > meanc + stddiv * stdc);
            abslst(i) = (ccva(loc(i)) > curve_threshold);
        end
        
        % =================================== %  
        % Approach 2 [outdated]
        % incircle and circumscribed circle
        % =================================== %  

        % mean square circle fit
        [xo,yo,R] = fr.circlefit(backpt(:,1),backpt(:,2));
        [xms,yms] = fr.mkcircle(xo,yo,R,bkt);

        % loader
        xmid = xms;
        ymid = yms;

        % calculation
        inlst = zeros(1,length(xmid));
        for i = 1:length(xmid)
            in = inpolygon(xmid(i),ymid(i),backpt(:,1),backpt(:,2));
            inlst(i) = in;
        end
        intseclst = zeros(1,length(inlst)-2);
        for i = 2:length(inlst)-1
            if inlst(i-1) ~= inlst(i) && inlst(i) == inlst(i+1)
                intseclst(i-1) = i;
            end
        end
        intseclst = intseclst(intseclst ~= 0);
        intseclst = [intseclst,intseclst(1)]; % close the intersections
    
        difflst = zeros(1,length(intseclst)-1);
        meandist = zeros(1,length(intseclst)-1);
        maxdist = zeros(1,length(intseclst)-1);
        stddist = zeros(1,length(intseclst)-1);
        areadist = zeros(length(intseclst)-1);
        diff2lst = zeros(length(intseclst)-1);
        for i = 1:length(intseclst)-1
            st = intseclst(i);
            et = intseclst(i+1);
            if et > st
                tdf = et - st;
                allptint = [xmid(st+1:et-1);ymid(st+1:et-1)]';
            else
                ll = length(xmid);
                tdf = et + (ll - st);
                firstpt = [xmid(st+1:ll);ymid(st+1:ll)]';
                secondpt = [xmid(1:et-1);ymid(1:et-1)]';
                allptint = [firstpt;secondpt];
            end
            distcnt = zeros(1,size(allptint,1));
            for j = 1:size(allptint,1)
                % calculate the distance and store in the set
                insidept = allptint(j,:); 
                [d,~,~] = fr.p_poly_dist(insidept(1),insidept(2),...
                    backpt(:,1),backpt(:,2));
                distcnt(j) = abs(d);
            end
            % store distance-related quantity
            difflst(i) = tdf; 
            meandist(i) = mean(distcnt);
            maxdist(i) = max(distcnt); 
            stddist(i) = std(distcnt);

            % calculate the embedded area (badly formulated)
            k1 = dsearchn(backpt,allptint(1,:));
            k2 = dsearchn(backpt,allptint(end,:));
            
            if k2-k1<0                
                sspt = [backpt(k1:end,:);backpt(1:k2,:)];
            else
                sspt = backpt(k1:k2,:);
            end
            
            ffpt = flip(allptint);
            formpt = [ffpt;sspt];
            formpt = [formpt;formpt(1,:)];

            closearea = polyarea(formpt(:,1),formpt(:,2));

            % store embedded area
            areadist(i) = closearea;
            diff2lst(i) = size(formpt,1) - size(allptint,1);    

        end
    
        intseclst(end) = [];
        selectedlst = intseclst;
    
        % some raw criteria
        for i = 1:length(selectedlst)
            if difflst(i) < 0.1 * bkt || maxdist(i) < 0.1 || meandist(i) < 0.04
                selectedlst(i) = 0;
            end
        end
    
        selectedlst(selectedlst == 0) = [];

        % plot
        adjustpeakcntt = 0;
        selectioncnt = 0;
    
        plt_crvind = 1;
        if plt_crvind == 1
            figure('units','normalized','outerposition',[0 0 1 1]);

            % plot reference figure of shape extraction
            subplot(1,4,1)
            imshow(orig);
            hold on
            plot(expixel,eypixel,'Color','r','LineWidth',1)
            hold off

            % plot reference figure for Approach 1
            subplot(1,4,2)
            plot(L2,ccva)
            hold on
            for i=1:length(peaks)
                if(indlst(i) == 1 && stdlst(i) == 1 && abslst(i) == 1)
                    scatter(L2(loc(i)),ccva(loc(i)),'red','filled')
                end
            end
            xlabel('s (arclength)')
            ylabel('\kappa (curvature)')
            title('Reference Plot for Approach 1')
    
            % plot Approach 1
            subplot(1,4,3)
            plot(backpt(:,1),backpt(:,2),'*')
            hold on
            plot(backex,-backey,'+')
            scatter(backpt(1,1),backpt(1,2),'black','filled')
            for i = 1:length(peaks)
                if(indlst(i) == 1 && stdlst(i) == 1 && abslst(i) == 1)
                    peakvalue = ccva(loc(i)); 
                    scatter(backpt(loc(i),1),backpt(loc(i),2),40,'cyan','filled')
                    adjustpeakcntt = adjustpeakcntt + 1;
                end
            end
            scatter(backcx,backcy,'red','filled')
            yline(-natpara.wcs * natpara.rr / 100)
            title(['Approach 1: Count = ', num2str(adjustpeakcntt)],...
                'FontSize',14)
            axis equal

            % plot Approach 2
            subplot(1,4,4)
            plot(backpt(:,1),backpt(:,2),'*')
            hold on
            plot(xmid,ymid,'LineStyle','--','LineWidth',1)
            for i = 1:length(intseclst)
                ii = intseclst(i);
                if sum(selectedlst == ii) == 0
                    scatter(xmid(ii),ymid(ii),50,'green','filled')
                end
            end
            scatter(xmid(1),ymid(1),70,'cyan','filled') % start point
            for i = 1:length(selectedlst)
                ii = selectedlst(i);
                scatter(xmid(ii),ymid(ii),50,'red','filled')
                selectioncnt = selectioncnt + 1;       
            end
            scatter(backcx,backcy,'filled')
            axis equal
            yline(-natpara.wcs * natpara.rr / 100)
            title(['Approach 2: Count = ', num2str(selectioncnt)],...
                'FontSize',14)

    
            sgtitle(["Second" string(f/rfr)]);
            close all;
        end
       
        % store the values
        adjustedcurvatureSet(iiii) = adjustpeakcntt;
        curvature3Set(iiii) = selectioncnt;
    
    end

    %% Polygon rotation (counterclockwise) [outdated]
    % 1: boundary integral
    % 2: body integral
    % 3: both present, and show the comparison, not for test application
    
    % parameter setting up
    makeupdata = 0;
    makeuprr = 0.05; % radius, only affect the amplitude
    perturb = 0; % whether to perturb the shape using harmonic function
    perturbrange = [0,180]; % certain angle range to add perturbance
    perturbmag = makeuprr / 10; % maximum magnitude of perturbance
    indinterval = 1; % grid of rotation angle
    numofside = 5; % n
    kappa_sval = 30; % kappa magnitude (unit not standardized)
    inward = 0; % kappa sign 

    % run
    if makeupdata == 1
        % create circle
        [x,y] = fr.mkcircle(0,0,makeuprr,natpara.numpt-1);
        x = x'; y = y';
        pt = [x,y];
        finalL = 2 * pi * makeuprr;
    elseif makeupdata == 2
        % create regular n-gon
        % including perturbation using harmonic functions or constant-curvature curve
        [x,y,origx,origy,LL] = fr.mkngon(0,0,makeuprr,natpara.numpt-1,intp,numofside,perturb,...
            perturbrange,perturbmag,kappa_sval,inward);

        pt = [x,y];
        finalL = LL;
    end
    
    % register interpolation result
    interpolation_cell{iiii} = [ex,-ey];

    rotationind = 0;
    if rotationind == 1
        crosscnt = bdyn.calbd_integral(pt,finalL,orig,ex,ey,cx,cy,...
            natpara,fr);
    elseif rotationind == 2
        [crosscnt2,t2,m2,s2,hw2] = bdyn.calbody_integral(orig,ashape,natpara,fr,...
            pt.*natpara.rr,makeupdata);

    elseif rotationind == 3
         [crosscnt1,equi1,t1,m1,s1,hw1,xc1,yc1] = bdyn.calbd_integral(pt,finalL,orig,ex,ey,cx,cy,...
                natpara,fr,indinterval);
         intergrid = 1;
         [crosscnt2,equi2,t2,m2,s2,hw2,xc2,yc2,oldxc2] = bdyn.calbody_integral(orig,ashape,natpara,...
             pt.*natpara.rr,makeupdata,indinterval,hw1,intergrid);

          % compare the result
         figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])
         subplot(1,3,1)     
         hold on
         scatter(pt(:,1),pt(:,2))
         axis equal
         yline(hw1,'r')
         yline(hw2,'g')
         scatter(xc2,yc2)
         legend('\alpha boundary','boundary-wl','body-wl')

         hold off

         subplot(1,3,2)
         hold on
         plot(t1,m1,'LineWidth',2)
         plot(t2,m2,'LineWidth',2)
         legend('Boundary Integral','Body Integral')
         xline(0,'LineStyle','--','LineWidth',1,'DisplayName','')
         yline(0,'LineStyle','--','LineWidth',1,'DisplayName','')
         title('Torque')
         hold off

         subplot(1,3,3)
         hold on
         plot(t1,s1,'LineWidth',2)
         plot(t2,s2,'LineWidth',2)
         legend('Boundary Integral','Body Integral')
         xline(0,'LineStyle','--','LineWidth',1,'DisplayName','')
         yline(0,'LineStyle','--','LineWidth',1,'DisplayName','')
         title('Moved Distance')
         hold off 

         % result showcase
         T = array2table([[xc1,yc1];[oldxc2,yc2];[xc2,yc2]],'VariableNames',...
             {'X-center','Y-center'},'RowName',{'Boundary','Old Body','Body'}); 
         
         % torque at steady
         m1_val0 = m1(abs(t1)==0);
         m2_val0 = m2(abs(t2)==0);

         roundheight = round(-natpara.wcs,5);

         H = array2table([[m1_val0,hw1,roundheight];[m2_val0,hw2,roundheight]],'VariableNames',...
             {'Torque when \theta=0','Water Height','Calibration Height'},...
             'RowNames',{'Boundary','Body'});

         P = array2table([crosscnt1;crosscnt2],'VariableNames',...
             {'Equilibrium Points'},'RowNames',{'Boundary','Body'});
        

         disp(T)
         disp(H)
         disp(P)

         close all

         Hbodyset(iiii) = hw2;
         cross2set(iiii) = crosscnt2;
   end
    
    % update the iteration parameter
    frameLab(iiii) = f;
    timeLab(iiii) = v.CurrentTime;    
end

% assure only two photos between flips are selected
indc = 0;
for tttt = 1:length(indInv{jj})
    indc = indc + indInv{jj}{tttt};
end
assert(indc == 2)

end

% output
savename = foldername(2:end-1);

% plot(timeLab,adjustedcurvatureSet,'-o')
% hold on
% plot(timeLab,curvature3Set,'-o')
% legend('curvature','circle')

% save(savename,'timeLab','adjustedcurvatureSet','curvature3Set')
% save(['meltrate-',savename],'alltInv','allmc','allnormalmc')

% plt_umrate
% plt_flipframe(alltInv,pttracker,indInv,closeoverind,wlineset(1))



