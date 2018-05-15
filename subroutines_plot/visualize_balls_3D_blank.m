% JN Kather, NCT Heidelberg, 2017
% compatible with model 2.0(TU/IM/MP):     yes
% compatible with 3D:                      no

function visualize_balls_3D_blank(mySystem,cnst,allSummaries)

currentSummary = allSummaries{end}; % current summary is last state
mysqueeze = @(varargin) (varargin); % squeeze function

    for i=1:numel(allSummaries)
        plotTimeline(i).TU_Num = allSummaries{i}.TU_Num;
        plotTimeline(i).LowAnti = allSummaries{i}.immuneEscape.lowAntigenicity;
        plotTimeline(i).LowAdju = allSummaries{i}.immuneEscape.lowAdjuvanticity;
        plotTimeline(i).losses = allSummaries{i}.LossFun;
    end
    
if numel(mySystem.TU.TUcells>0) && numel(mySystem.IM.IMcells>0) ...
        && numel(mySystem.MP.MPcells>0) % check that all major players exist
    
    % prepare tumor antigenicity grid before cutting out tumor cell slice
    antigenicityGrid = double(zeros(size(mySystem.grid.L)));   % preallocate antigenicity grid
    antigenicityGrid(mySystem.TU.TUcells) = mySystem.TU.TUprop.Antigen;   % preallocate antigenicity grid
      
    figure(1)    % prepare figure 1 = tumor 3D plot
    clf('reset')
    % show full domain and add decorations
    xsize = 4; ysize = 10;
    posi.MainPanel = [1:6,11:16,21:26,31:36];
    posi.AntiPanel = 17:18;
    posi.AdjuPanel = 27;
    posi.AdjuPanel2 = 28;
    posi.PopPanel1 = 37;
    posi.PopPanel2 = 38;
    posi.TimelinePanel = 7:8;
    posi.lymTopo = 9:10; 
    posi.myeTopo = 19:20;
    posi.AntiAdjuPanel = 29:30;
    posi.lossPanel = 39:40;
    
    subplot(xsize,ysize,posi.MainPanel);
    
    TUcolors = double(hot(round(2*cnst.tumorColorLevels))); %colormap for tumor cells
    TUcolors = TUcolors(0.5*cnst.tumorColorLevels:1.5*cnst.tumorColorLevels,:);% crop colormap
    IMcolors = flipud(double(blugr(double(mySystem.params.IMkmax)+3))); % color map for lymphocytes
    MPcolors = [0 0.3 0; 0 0.8 0]; % macrophage color map, only two states
    
    N1 = mySystem.grid.N;
    N2 = mySystem.grid.M;
    N3 = mySystem.grid.P;
    
    %clearing cells from the boundary
    L = false(size(mySystem.grid.L));
    Lymph = L;
    Mac = L;
    Nec = L;
    Fib = L;
    
    clear90 = true;
    
    if clear90
    %clearing cells within 90 degree angle to have cutout
    [x, y, ~] = ind2sub([N1 N2 N3],mySystem.TU.TUcells);
    
    %calculating center of mass
    X = mean(x); Y = mean(y); %Z = mean(z);   
    indx = (x > X) & (y > Y);
    mySystem.TU.TUcells(indx) = []; % remove tumor cells
    mySystem.TU.TUprop.Antigen(indx) = []; % also remove antigenicity property
    
    [x, y, ~] = ind2sub([N1 N2 N3],mySystem.IM.IMcells);
    mySystem.IM.IMcells((x > X) & (y > Y)) = [];
    
    [x, y, ~] = ind2sub([N1 N2 N3],mySystem.MP.MPcells);
    mySystem.MP.MPcells((x > X) & (y > Y)) = [];

    necro = int32(find(mySystem.grid.Ln)');
    [x, y, ~] = ind2sub([N1 N2 N3],necro);
    necro((x > X) & (y > Y)) = [];
    
    fibr = int32(find(mySystem.grid.Lf)');
    [x, y, ~] = ind2sub([N1 N2 N3],fibr);
    fibr((x > X) & (y > Y)) = [];
    end
    
    % update grids
    L( mySystem.TU.TUcells) = true;
    Lymph( mySystem.IM.IMcells) = true;
    Mac(mySystem.MP.MPcells) = true;
    Nec(necro) = true;
    Fib(fibr) = true;

    
    blLevels = 1;
    if blLevels %if perform blur
        
        %auxilary variable with indices to the cell neighborhood
        aux = int32([[-N1-1 -N1 -N1+1 -1 1 N1-1 N1 N1+1] ...
            [-N1-1 -N1 -N1+1 -1 1 N1-1 N1 N1+1]-N1*N2 ...
            [-N1-1 -N1 -N1+1 -1 1 N1-1 N1 N1+1]+N1*N2])';
        
        S = [mySystem.TU.TUcells unique(reshape(bsxfun(@plus,mySystem.TU.TUcells,aux),1,[]))];
        S2 = bsxfun(@plus,S,aux);
        S2(S2<1) = []; S2(S2>N1*N2*N3) = [];
        
        SL = [mySystem.IM.IMcells unique(reshape(bsxfun(@plus,mySystem.IM.IMcells,aux),1,[]))];
        SL2 = bsxfun(@plus,SL,aux);
        SL2(SL2<1) = []; SL2(SL2>N1*N2*N3) = [];
        
        SM = [mySystem.MP.MPcells unique(reshape(bsxfun(@plus, mySystem.MP.MPcells,aux),1,[]))];
        SM2 = bsxfun(@plus,SM,aux);
        SM2(SM2<1) = []; SM2(SM2>N1*N2*N3) = [];
        
        SN = [necro unique(reshape(bsxfun(@plus, necro,aux),1,[]))];
        SN2 = bsxfun(@plus,SN,aux);
        SN2(SM2<1) = []; SN2(SN2>N1*N2*N3) = [];
        
        SF = [fibr unique(reshape(bsxfun(@plus, fibr,aux),1,[]))];
        SF2 = bsxfun(@plus,SF,aux);
        SF2(SF2<1) = []; SF2(SF2>N1*N2*N3) = [];
        
        %changing lattice from logical variable to float
        L = single(L);
        Lymph = single(Lymph);
        Mac = single(Mac);
        Nec = single(Nec);
        Fib = single(Fib);
        for i = 1:blLevels %for number of blurs
            L(S) = mean(L(S2)); %taking the average of neighborhood
            Lymph(SL) = mean(Lymph(SL2)); %taking the average of neighborhood
            Mac(SM) = mean(Mac(SM2)); %taking the average of neighborhood
            Nec(SN) = mean(Nec(SN2)); %taking the average of neighborhood
            Fib(SF) = mean(Fib(SF2)); %taking the average of neighborhood
        end
    end
    
    
    %calculating isosurfaces and plotting
    hold on
    %necrosis
    pNf = patch(isosurface(1:N1,1:N2,1:N3,Nec,0.5));
    isonormals(1:N1,1:N2,1:N3,Nec,pNf)
    set(pNf,'FaceColor','k','EdgeColor','none');
    
    %fibrosis
    pLf = patch(isosurface(1:N1,1:N2,1:N3,Fib,0.5));
    isonormals(1:N1,1:N2,1:N3,Fib,pLf)
    set(pLf,'FaceColor',[128 128 128]/255,'EdgeColor','none'); % 128 128 128
    
    % lymphocytes
    pL = patch(isosurface(1:N1,1:N2,1:N3,Lymph,0.25));
    isonormals(1:N1,1:N2,1:N3,Lymph,pL)
    set(pL,'FaceColor',IMcolors(end,:),'EdgeColor','none');
    
    % macrophages
    pM = patch(isosurface(1:N1,1:N2,1:N3,Mac,0.25));
    isonormals(1:N1,1:N2,1:N3,Mac,pM)
    set(pM,'FaceColor',MPcolors(1,:),'EdgeColor','none');
    
    
    %[faces,verts,colors] = isosurface(1:N1,1:N2,1:N3,L,0.25,antigenicityGrid); % show tumor surface in antigen colors
  %     patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
%     'FaceColor',TUcolors(1,:),'EdgeColor','none');
%     
    
     % plot tumor as one color only
     pT = patch(isosurface(1:N1,1:N2,1:N3,L,0.25));
     isonormals(1:N1,1:N2,1:N3,L,pT)
     isonormals(1:N1,1:N2,1:N3,L,pT)
     set(pT,'FaceColor',TUcolors(1,:),'EdgeColor','none');
    
    hold off
    % add time counter to image
    title([num2str(mySystem.grid.StepsDone/2),' days'],...
        'Color','k','FontWeight','bold','FontSize',18,'VerticalAlignment','top')
    set(gcf,'Color','w');
    axis equal

    xlim([1 N1]);
    ylim([1 N2]);
    zlim([1 N3]);
   
    view([66 43]); %135 20
    camlight
    camlight('right') % add more light
    camlight('left')  % add more light
    lighting gouraud
    
    % prepare decorations
    myPos = get(0,'Screensize');
    myPos(1) = myPos(1) + myPos(3)*(1/10);
    myPos(2) = myPos(2) + myPos(4)*(1/10);
    myPos(3) = myPos(3)*(8/10);
    myPos(4) = myPos(4)*(8/10);
    set(gcf, 'Position', myPos);  % enlarge figure
     
    % antigenicity panel
    subplot(xsize,ysize,posi.AntiPanel);
    if(numel(mySystem.TU.TUprop.Antigen)>3)
        h = histogram(mySystem.TU.TUprop.Antigen,0:0.1:1);
        h.EdgeColor = 'none';
        h.FaceColor = 'k';
        %plot(ky,kx/max(kx),'k','LineWidth',1.5);
    end
    hold on
    currax = axis();
    plot([mySystem.params.antiThresh,mySystem.params.antiThresh],[0,currax(4)],'k:','LineWidth',1.5);
    %axis([0,cnst.maxAntigenicity,0,1])
    title('antigenicity');
    
    % adjuvanticity panel
    subplot(xsize,ysize,posi.AdjuPanel);
    if(numel(mySystem.TU.TUcells)>3)
        [x, y, z] = ind2sub([mySystem.grid.N mySystem.grid.M mySystem.grid.P],mySystem.TU.TUcells);
        %calculating center of mass
        X = mean(x); Y = mean(y); Z = mean(z);
        pS = slice(1:mySystem.grid.N,1:mySystem.grid.M,1:mySystem.grid.P,(mySystem.grid.AdjuMap+10)/20,X,Y,Z);
        set(pS,'EdgeColor','none');
        xlim([1 mySystem.grid.N]);
        ylim([1 mySystem.grid.M]);
        zlim([1 mySystem.grid.P]);
    end
    view([66 43]);
    caxis([0 1]);
    title('adjuvanticity');
    
    % antigenicity 3D panel
    subplot(xsize,ysize,posi.AdjuPanel2);
    if(numel(mySystem.TU.TUcells)>3)
        [x, y, z] = ind2sub([mySystem.grid.N mySystem.grid.M mySystem.grid.P],mySystem.TU.TUcells);
        %calculating center of mass
        X = mean(x); Y = mean(y); Z = mean(z);
        pS = slice(1:mySystem.grid.N,1:mySystem.grid.M,1:mySystem.grid.P,antigenicityGrid,X,Y,Z);
        set(pS,'EdgeColor','none');
        xlim([1 mySystem.grid.N]);
        ylim([1 mySystem.grid.M]);
        zlim([1 mySystem.grid.P]);
    end
    colormap hot
    %colorbar;
    caxis([0 1]);
    view([66 43]);
    title('antigenicity');
    

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
    plot([1 0.5*(cnst.drawWhen:cnst.drawWhen:(cnst.drawWhen*numel(plotTimeline)))],...
        [1 cell2mat(mysqueeze(plotTimeline.TU_Num))],'k','LineWidth',1.5);
    end
    axis([0 0.5*cnst.nSteps 0 numel(mySystem.grid.L)/2])
    ylabel('tumor cells');
    title('tumor timeline');
    
    % AntiAdjuPanel
    subplot(xsize,ysize,posi.AntiAdjuPanel)
    if ~isempty(plotTimeline)
    xanti = cell2mat(mysqueeze(plotTimeline.LowAnti));
    yadju = cell2mat(mysqueeze(plotTimeline.LowAdju));
    plot(xanti, yadju,'k','LineWidth',1.1);
    hold on
    scatter(xanti(1),yadju(1),15,'ro','filled'); % label the beginning
    hold off
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
    axis([0 0.5*cnst.nSteps 0 5]);
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
    
   % print(gcf,[num2str(round(rand()*100000)),'_snap.png'],'-dpng','-r900');
    % show parameters
    paramtext1 = ['TUps = ',num2str(mySystem.params.TUps),...
                 ' | TUpmig = ',num2str(mySystem.params.TUpmig),...
                 ' | TUpmut = ',num2str(mySystem.params.TUpmut),...
                 ' | initialSeed = ',num2str(mySystem.params.initialSeed),...
                 ' | probSeedFibr = ',num2str(mySystem.params.probSeedFibr),...
                 ...
                 ' | IMinfluxProb = ',num2str(mySystem.params.IMinfluxProb),...
                 ' | IMinfluxRate = ',num2str(mySystem.params.IMinfluxRate),...
                 ' | IMrateDynamic = ',num2str(mySystem.params.IMrateDynamic),...
                 ' | IMinfluxModifier = ',num2str(mySystem.params.IMinfluxModifier),...
                 ' | IMpkill = ',num2str(mySystem.params.IMpkill),...
                 10,...
                 ' | MPinfluxProb = ',num2str(mySystem.params.MPinfluxProb),...
                 ' | MPinfluxRate = ',num2str(mySystem.params.MPinfluxRate),...
                 ' | MPrateDynamic = ',num2str(mySystem.params.MPrateDynamic),...
                 ' | MPinfluxModifier = ',num2str(mySystem.params.MPinfluxModifier),...
                 ...
                 ' | stromaPerm = ',num2str(mySystem.params.stromaPerm),...
                 ' | p-loss = ',num2str(currentSummary.p_losses),...
                 ' | m-loss = ',num2str(currentSummary.LossFun)];
             
    s = suptitle(paramtext1);
    s.FontSize = 10;
    
else
    disp('no tumor cells, could not plot anything');
end

end