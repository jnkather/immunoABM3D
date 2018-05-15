% JN Kather 2017, jakob.kather@nct-heidelberg.de
%
% compatible with model 2.0(TU/IM/MP):     yes
% compatible with 3D:                      no

function summaryOut = summarizeSystem_2D(mySystem,cnst)
%summarizeSystem_NM summarizes the state of a 2D system
%   Input is the "mySystem" structure with fields grid, params, TU, IM
%   and the "cnst" structure, containing basic constants
%   This function genreates a structure "summaryOut" containing all 
%   relevant measurements of the system's state

    % create basic results
    summaryOut.TU_Num = numel(mySystem.TU.TUcells); % tumor cell number
    summaryOut.TU_FracStem = sum(mySystem.TU.TUprop.isStem) ...
        / numel(mySystem.TU.TUprop.isStem); % stem cell fraction
    summaryOut.IM_Num = numel(mySystem.IM.IMcells); % immune cell number
    summaryOut.IM_FracExhaust = sum(mySystem.IM.IMprop.Kcap == 0) ...
        / numel(mySystem.IM.IMcells); % fraction of exhausted immune cells
    
	% if there is a chemotaxis map, then create spatial results
    if isfield(mySystem.grid,'ChtaxMap')
    ChtaxMap2 = double(bwdist(mySystem.grid.Lt,'euclidean')); 
    Mask_fullTumorAux = imfill(ChtaxMap2<cnst.inTumor,'holes');   % binary mask ROI 1
    Mask_tumorCore = imerode(Mask_fullTumorAux,strel('disk',round(cnst.marginSize/2),6));
    Mask_marginIn = Mask_fullTumorAux & ~Mask_tumorCore;
    Mask_marginOut = ChtaxMap2<cnst.marginSize & ~Mask_fullTumorAux; % binary mask ROI 2
    Mask_distantOut = ChtaxMap2<cnst.around & ~Mask_marginOut & ~Mask_fullTumorAux;  % binary mask ROI 3
	
	% count lymphocyte (IM) cells in regions PER GRID CELL
    summaryOut.IM_tumorCore = sum(Mask_tumorCore(mySystem.IM.IMcells))/sum(Mask_tumorCore(:));
    summaryOut.IM_marginIn = sum(Mask_marginIn(mySystem.IM.IMcells))/sum(Mask_marginIn(:));
    summaryOut.IM_marginOut = sum(Mask_marginOut(mySystem.IM.IMcells))/sum(Mask_marginOut(:));
    summaryOut.IM_distantOut = sum(Mask_distantOut(mySystem.IM.IMcells))/sum(Mask_distantOut(:));
	
    % count macrophage (MP) cells in regions PER GRID CELL
    summaryOut.MP_tumorCore = sum(Mask_tumorCore(mySystem.MP.MPcells))/sum(Mask_tumorCore(:));
    summaryOut.MP_marginIn = sum(Mask_marginIn(mySystem.MP.MPcells))/sum(Mask_marginIn(:));
    summaryOut.MP_marginOut = sum(Mask_marginOut(mySystem.MP.MPcells))/sum(Mask_marginOut(:));
    summaryOut.MP_distantOut = sum(Mask_distantOut(mySystem.MP.MPcells))/sum(Mask_distantOut(:));
    
	% more features: tumor/stroma ratio, tumor/necrosis ratio, stromal lymphocytes fraction
    summaryOut.TU_Stro_ratio_log = log(double(numel(mySystem.TU.TUcells))/double(sum(mySystem.grid.Lf(:))));
    summaryOut.Stro_Fraction = double(sum(mySystem.grid.Lf(:)))/double(numel(mySystem.TU.TUcells));
    summaryOut.TU_Necr_ratio_log = log(double(numel(mySystem.TU.TUcells))/double(sum(mySystem.grid.Ln(:))));
    summaryOut.IM_instroma = sum(mySystem.grid.Lf(mySystem.IM.IMcells))/numel(mySystem.IM.IMcells);
    
    % finally: calculate the formal loss function
    [summaryOut.p_losses, summaryOut.LossFun] = lossFunction(mySystem,cnst);
    
    % calculate immune escape indices: fraction of tumor cells in
    % low-antigenicity (sub-threshold) or low-adjuvanticity (sub-threshold)
    % conditions
    adjuvanticities =  mySystem.grid.AdjuMap(mySystem.TU.TUcells); % get all tumor cell adjuvanticities
    summaryOut.immuneEscape.lowAntigenicity  = sum(mySystem.TU.TUprop.Antigen <= mySystem.params.antiThresh) / numel(mySystem.TU.TUcells);
    summaryOut.immuneEscape.meanAntigenicity = mean(mySystem.TU.TUprop.Antigen);
    summaryOut.immuneEscape.lowAdjuvanticity = sum( adjuvanticities <= mySystem.params.adjuThresh)  / numel(mySystem.TU.TUcells);
    
    disp(['immune escape measurements: ',10,'low anti: ',num2str(round(100*summaryOut.immuneEscape.lowAntigenicity)),'%',...
        ', low adju: ', num2str(round(100*summaryOut.immuneEscape.lowAdjuvanticity)),'%']);
    end
    
    % copy hyper-parameters
    summaryOut.stepsDone = mySystem.grid.StepsDone;
    summaryOut.N = mySystem.grid.N;
    summaryOut.M = mySystem.grid.M;
    
end