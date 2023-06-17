close all

dataregs = cell(1,length(allmc));

for j=1:length(allmc)
    mc = allmc{j};
    time = alltInv{j}/rfr;
    centers = centerlst{j};
    timediff = time(end)-time(1); % s
    normalmc = allnormalmc{j};
    % initial shape
    startpt = mc{1};
    startcenter = centers{1};
    % normal vectors subject to the initial shape
    startnormal = normalmc{1};
    % second shape
    endpt = mc{end};
    endcenter = centers{end};

    stx = startpt(:,1); sty = startpt(:,2); % startpt [fixed]
    [~,deeppt] = min(sty);
    bottomx = stx(deeppt);

    colorlst = ["g","b","c","m","#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30"];

    % manipulate the frames
    relshift = endcenter(:,2) - startcenter(:,2);
    coffstd = 3;
    shiftlst = linspace(-coffstd*relshift,coffstd*relshift,coffstd*2+1);

    shitarcrate = cell(1,length(shiftlst));
    shiftgraph = cell(1,length(shiftlst)+1);

    for iN = 1:length(shiftlst)
        legendCell{iN} = num2str(shiftlst(iN),"Move=%-d m");
    end

    figure('units','normalized','outerposition',[0 0 1 1])
    tiledlayout(2,1)
    ax1 = nexttile;
    xlabel(ax1,'Cumulative Arclength')
    ylabel(ax1,'Melting Rate')
    title(ax1,['Between ', num2str(j), ' and ', num2str(j+1),',center of mass difference=', ...
        num2str(relshift),'m'])
    hold(ax1,'on')

    ax2 = nexttile;
    
    plot(ax2,stx,sty,'o','Color','r');
    % save first shape
    shiftgraph{1} = [stx,sty];

    title(ax2,'Graph Illustration')
    axis equal
    hold(ax2,'on')
   
    for kkk = 1:length(shiftlst)
        testmove = shiftlst(kkk);
        disp(['Move=',num2str(testmove)])
        clear newendpt
        newendpt(:,2) = endpt(:,2) + testmove; % move y
        newendpt(:,1) = endpt(:,1); % keep x

        plot(ax2,newendpt(:,1),newendpt(:,2),'*','Color',colorlst(kkk))
        % save movement shape
        shiftgraph{kkk+1} = newendpt;
    
        N = length(stx);
        Kgr = 1501;
    
        trailrange = linspace(-5*10^(-5),10*10^(-5),Kgr);
        meltratereg = zeros(1,N);
        parfor (i = 1:N,10)
            for k = 1:Kgr
                meltrate = trailrange(k);
                newhypt = [stx(i)+meltrate*timediff*startnormal(i,1);...
                    sty(i)+meltrate*timediff*startnormal(i,2)];
                [in,on] = inpolygon(newhypt(1),newhypt(2),newendpt(:,1),newendpt(:,2));
                if in == 1
                    meltratereg(i) = meltrate;
                    break
                end
            end
        end
    
        ds = sqrt((stx(2:end)-stx(1:end-1)).^2+(sty(2:end)-sty(1:end-1)).^2);
        s = zeros(N,1); % arclength accumulation
        for n = 1:N-1
            s(n+1) = s(n) + ds(n);
        end
        finalL = s(end);
        
        plot(ax1,s,meltratereg,'o','Color',colorlst(kkk));
        % save movement analysis result
        shitarcrate{kkk} = [s,meltratereg'];
    end
    legend(ax1,legendCell);
    saveas(gcf,['move-vertical-',num2str(j),'.jpg'])

    save('arc_rate.mat','shitarcrate');
    save('mov_graph.mat','shiftgraph');
    close all

    dataregs{j}(:,1) = s/finalL;
    dataregs{j}(:,2) = meltratereg;
    disp([num2str(j),' flip is calculated.'])

end

figure()
hold on
xxx = linspace(0,1,100000);
kk = zeros(1,length(xxx));
for i = 1:length(dataregs)
    plot(dataregs{i}(:,1),dataregs{i}(:,2),'-');
    vq1 = interp1(dataregs{i}(:,1),dataregs{i}(:,2),xxx);
    kk = kk + vq1;
end
kk = kk/length(dataregs);
plot(xxx,kk,'-o','LineWidth',2);
xlabel('Proportion of Cumulative Arclength','FontSize',14)
ylabel('Melting Rate (m/s)','FontSize',14)
saveas(gcf,['meltrate-',foldername(1:end-1),'-end','.jpg'])



