% JN Kather, NCT Heidelberg, 2017
% compatible with model 2.0(TU/IM/MP):     yes
% compatible with 3D:                      no

function visualize_balls_2D_blank(mySystem,cnst,allSummaries)

currentSummary = allSummaries{end}; % current summary is last state
mysqueeze = @(varargin) (varargin); % squeeze function

    for i=1:numel(allSummaries)
        plotTimeline(i).TU_Num = allSummaries{i}.TU_Num;
        plotTimeline(i).LowAnti = allSummaries{i}.immuneEscape.lowAntigenicity;
        plotTimeline(i).LowAdju = allSummaries{i}.immuneEscape.lowAdjuvanticity;
        plotTimeline(i).losses = allSummaries{i}.LossFun;
    end

    if numel(mySystem.TU.TUcells>0)
        
    % prepare figure
    clf('reset')
    % show full domain and add decorations
    xsize = 4; ysize = 10;
    posi.MainPanel = [1:6,11:16,21:26,31:36];
    posi.AntiPanel = 17:18;
    posi.AdjuPanel = 27:28;
    posi.PopPanel1 = 37;
    posi.PopPanel2 = 38;
    posi.TimelinePanel = 7:8;
    posi.lymTopo = 9:10; 
    posi.myeTopo = 19:20;
    posi.AntiAdjuPanel = 29:30;
    posi.lossPanel = 39:40;
    subplot(xsize,ysize,posi.MainPanel);
    
    % create color vectors for all cell types for all cell states
    % % VISUALIZE REMAINING PROLIFERATION CYCLES OR ANTIGENICITY
    TUcolors = double(hot(round(2*cnst.tumorColorLevels))); %colormap for tumor cells 
    TUcolors = TUcolors(0.5*cnst.tumorColorLevels:1.5*cnst.tumorColorLevels,:);% crop colormap
    IMcolors = flipud(double(blugr(double(mySystem.params.IMkmax)+3))); % color map for lymphocytes
    MPcolors = [0 0.3 0; 0 0.8 0]; % macrophage color map, only two states
    
    % prepare figure background
    backgroundImage = 255* double(ones(mySystem.grid.N,mySystem.grid.M,3));
    ch1 = backgroundImage(:,:,1);
    ch2 = backgroundImage(:,:,2);
    ch3 = backgroundImage(:,:,3);
    
    % add STROMA to background image
    ch1(mySystem.grid.Lf) = 128/255;
    ch2(mySystem.grid.Lf) = 128/255;
    ch3(mySystem.grid.Lf) = 128/255;
    
    % add necrosis to background image
    ch1(mySystem.grid.Ln) = 0;
    ch2(mySystem.grid.Ln) = 0;
    ch3(mySystem.grid.Ln) = 0;
    
    % add 2 mm scale bar to background image
    beginSc = 10;
    width = 3; 
    len = 134; % 2 mm length
    ch1((end-beginSc-width):(end-beginSc),beginSc:(beginSc+len)) = 0;
    ch2((end-beginSc-width):(end-beginSc),beginSc:(beginSc+len)) = 0;
    ch3((end-beginSc-width):(end-beginSc),beginSc:(beginSc+len)) = 0;
    
    % collect background image channels
    backgroundImage(:,:,1) = ch1;
    backgroundImage(:,:,2) = ch2;
    backgroundImage(:,:,3) = ch3;
                                     
    imshow(backgroundImage); % show background
    hold on % prepare to plot on top
   
    % prepare decorations
    myPos = get(0,'Screensize');
    myPos(1) = myPos(1) + myPos(3)*(1/10);
    myPos(2) = myPos(2) + myPos(4)*(1/10);
    myPos(3) = myPos(3)*(8/10);
    myPos(4) = myPos(4)*(8/10);
    set(gcf, 'Position', myPos);  % enlarge figure
    
    % create tumor cell coordinates
    ytu = mod(mySystem.TU.TUcells,mySystem.grid.N);
    xtu = mySystem.TU.TUcells/(mySystem.grid.N-1);   
    % retrieve tumor cell colors
    % myTUc = TUcolors(mySystem.TU.TUprop.Pcap+5,:); % REMAINING PROLIF
    try % try to draw tumor, skip if this is impossible
        myTUc = TUcolors(1+round(mySystem.TU.TUprop.Antigen/cnst.maxAntigenicity*cnst.tumorColorLevels),:); % ANTIGENICITY
        % plot tumor cells
        scatter(xtu,ytu,mySystem.params.mycellsize,myTUc,'filled')
    catch
        max(mySystem.TU.TUprop.Antigen+1)
        whos
        error('failed to draw tumor cells');       
    end
       
    try % try to draw lymphocytes, skip if this is impossible
    % create lymphocyte coordinates
    yim = mod(mySystem.IM.IMcells,mySystem.grid.N);
    xim = mySystem.IM.IMcells/(mySystem.grid.N-1);
    % retrieve lymphocyte colors
    myIMc = IMcolors(mySystem.IM.IMprop.Kcap+1,:);
    % plot lymphocytes
    scatter(xim,yim,mySystem.params.mycellsize,myIMc,'filled')
    catch
        disp('no lymphocytes could be plotted');
    end
    
    try % try to draw macrophages, skip if this is impossible
    % create macrophages coordinates
    ymp = mod(mySystem.MP.MPcells,mySystem.grid.N);
    xmp = mySystem.MP.MPcells/(mySystem.grid.N-1);
    % retrieve macrophages colors
    myMPc = MPcolors(mySystem.MP.MPprop.State,:);
    % plot macrophages
    scatter(xmp,ymp,mySystem.params.mycellsize,myMPc,'filled')
    catch
        disp('no macrophages could be plotted');
    end 
    
    % show full domain and add decorations
    set(gca,'XLim',[0 mySystem.grid.M]);
    set(gca,'YLim',[0 mySystem.grid.N]);
    axis equal off
    set(gcf,'Color','w');
    
    % add time counter to image
    text(beginSc,beginSc/2,[num2str(mySystem.grid.StepsDone/2),' days'],...
        'Color','k','FontWeight','bold','FontSize',18,'VerticalAlignment','top')
    
    % show parameters
    paramtext1 = ['TUps = ',num2str(mySystem.params.TUps),...
                 ' | TUpmig = ',num2str(mySystem.params.TUpmig),...
                 ' | TUpmut = ',num2str(mySystem.params.TUpmut),...
                 ' | TUdadju = ',num2str(mySystem.params.TUdadju),...
                 ' | initialSeed = ',num2str(mySystem.params.initialSeed),...
                 ' | probSeedFibr = ',num2str(mySystem.params.probSeedFibr),...
                 ' | necThresh = ',num2str(mySystem.params.necThresh),10,...
                 ' IMinfluxProb = ',num2str(mySystem.params.IMinfluxProb),...
                 ' | IMinfluxRate = ',num2str(mySystem.params.IMinfluxRate),...
                 ' | IMinfluxModifier = ',num2str(mySystem.params.IMinfluxModifier),...
                 ' | IMpkill = ',num2str(mySystem.params.IMpkill),...
                 ' | MPinfluxProb = ',num2str(mySystem.params.MPinfluxProb),...
                 ' | MPinfluxModifier = ',num2str(mySystem.params.MPinfluxModifier),10,...
                 'MPprepola = ',num2str(mySystem.params.MPprepola),...
                 ' | adjuDecay = ',num2str(mySystem.params.adjuDecay),...
                 ' | stromaPerm = ',num2str(mySystem.params.stromaPerm),10,...
                 'loss = ',num2str(currentSummary.LossFun),...
                 ' | lowAnti% = ',num2str(round(100*currentSummary.immuneEscape.lowAntigenicity)),...
                 ' | lowAdju% = ',num2str(round(100*currentSummary.immuneEscape.lowAdjuvanticity))];
    text(beginSc-30,-50,paramtext1,'FontSize',11);
    
    hold off
    
    % add side panels
    
    % antigenicity panel
    subplot(xsize,ysize,posi.AntiPanel); 
    if(numel(mySystem.TU.TUprop.Antigen)>3)    
    [kx,ky] = ksdensity(mySystem.TU.TUprop.Antigen,'width',cnst.antigenKernelWidth); 
    plot(ky,kx/max(kx),'k','LineWidth',1.5);
    end
    hold on
    plot([mySystem.params.antiThresh,mySystem.params.antiThresh],[0,1],'k:','LineWidth',1.5);
    axis([0,cnst.maxAntigenicity,0,1])
    title('antigenicity');
    
    % adjuvanticity panel
    subplot(xsize,ysize,posi.AdjuPanel);
    if(numel(mySystem.TU.TUcells)>3)  
    imagesc(mySystem.grid.AdjuMap),axis equal tight off,
    end
    colormap(gray),colorbar,
    caxis([-10 10]);
    title('adjuvanticity');
    
    % population panel 1/2
    subplot(xsize,ysize,posi.PopPanel1)
    popcount1 = [sum(mySystem.TU.TUprop.isStem),...
                sum(~mySystem.TU.TUprop.isStem),...
                sum(mySystem.grid.Lf(:)),...
                sum(mySystem.grid.Ln(:))];
    poplabel1 = categorical({'TU_{stem}','TU_{norm}','STRO','NEC'});
    bar(poplabel1,popcount1);
    set(gca,'XTickLabelRotation',90);
    
    % population panel 2/2
    subplot(xsize,ysize,posi.PopPanel2)
    popcount2 = [sum(mySystem.IM.IMprop.Kcap>0),...
                sum(mySystem.IM.IMprop.Kcap==0),...
                sum(mySystem.MP.MPprop.State==1),...
                sum(mySystem.MP.MPprop.State==2)];
    poplabel2 = categorical({'LYM_{act}','LYM_{exh}','MP_{naive}','MP_{pro}'});
    bar(poplabel2,popcount2);
    set(gca,'XTickLabelRotation',90);
    
    % TU Num timeline Panel
    subplot(xsize,ysize,posi.TimelinePanel)
    if ~isempty(plotTimeline)
    plot([0 0.5*(cnst.drawWhen:cnst.drawWhen:(cnst.drawWhen*numel(plotTimeline)))],...
        [1 cell2mat(mysqueeze(plotTimeline.TU_Num))],'k','LineWidth',1.5);
    end
    axis([0 0.5*cnst.nSteps 0 numel(mySystem.grid.L)/5])
    ylabel('tumor cells');
    title('tumor timeline');
    
    % AntiAdjuPanel
    subplot(xsize,ysize,posi.AntiAdjuPanel)
    if ~isempty(plotTimeline)
    plot(cell2mat(mysqueeze(plotTimeline.LowAnti)),...
        cell2mat(mysqueeze(plotTimeline.LowAdju)),'k','LineWidth',1.1);
    end
    xlabel('low anti %'); ylabel('low adju %');
    axis equal
    axis([0 1 0 1]);
    title('immune escape');
    
    % Loss Panel
    subplot(xsize,ysize,posi.lossPanel)
    if ~isempty(plotTimeline)
    plot([0 0.5*(cnst.drawWhen:cnst.drawWhen:(cnst.drawWhen*numel(plotTimeline)))],...
        [0 cell2mat(mysqueeze(plotTimeline.losses))],'k','LineWidth',1.5);
    end
    xlabel('time (days)'); ylabel('mean loss');
    axis square
    axis([0 0.5*cnst.nSteps 0 10]);
    title('loss'); 
    

    % normal color for topography
    colornorm=cnst.topoColorNorm;
    xbase = 1;
    ybase = 1;
    
    % lymphocytes topography panel
    subplot(xsize,ysize,posi.lymTopo)
    colordata.TU_CORE = currentSummary.IM_tumorCore;
    colordata.MARG_500_IN = currentSummary.IM_marginIn;
    colordata.MARG_500_OUT = currentSummary.IM_marginOut;
    mytitle = 'lymphoid topography';
    drawcircles(xbase,ybase,colordata,mytitle,colornorm)
    axis square off
    %title('lymphoid topography');
    
    % macrophages topography panel
    subplot(xsize,ysize,posi.myeTopo)
    colordata.TU_CORE = currentSummary.MP_tumorCore;
    colordata.MARG_500_IN = currentSummary.MP_marginIn;
    colordata.MARG_500_OUT = currentSummary.MP_marginOut;
    mytitle = 'myeloid topography';
    drawcircles(xbase,ybase,colordata,mytitle,colornorm)
    axis square   off
    %title('myeloid topography');
    
    else
        disp('no tumor cells, could not plot anything');
    end

end