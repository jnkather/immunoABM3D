% JN Kather 2018
% this function calculates the loss of the system, i.e. the deviation from
% a desired (realistic) state. This desired state is defined by actual
% clinical data, for example TCGA.

function [p_loss, m_loss] = lossFunction(mySystem,cnst)
% return partial losses and mean loss
if isfield(cnst,'lossFunction')
    switch cnst.lossFunction
        case 'default' % default loss function does not contain tumor stem cell fraction. Should be used for 2D.
            
        % -- define desired values
        if cnst.verbose
            disp('calculating default loss function');
        end
%         % stroma-tumor ratio, from TCGA, N=236 pts.
%         stroma_tumor_ratio_mean = 0.1992;
%         stroma_tumor_ratio_std  = 0.4547;
%         % lymphocyte-tumor ratio, own data, N=15 pts.
%         lympho_tumor_ratio_mean = 3.72/100;
%         lympho_tumor_ratio_std  = 3.13/100;
%         % macrophage-tumor ratio, own data, N=15 pts.
%         macro_tumor_ratio_mean  = 17.04/100;
%         macro_tumor_ratio_std   = 12.92/100;

        % stroma-tumor ratio, from own data, N=14 patients, 16 Feb 2018
        stroma_tumor_ratio_mean = 0.2263;
        stroma_tumor_ratio_std  = 0.1562;
        % lymphocyte-tumor ratio, from own data, N=14 patients, 16 Feb 2018
        lympho_tumor_ratio_mean = 0.0332;
        lympho_tumor_ratio_std  = 0.0411;
        % macrophage-tumor ratio, from own data, N=14 patients, 16 Feb 2018
        macro_tumor_ratio_mean  = 0.0978;
        macro_tumor_ratio_std   = 0.0825;

        % exhausted-total lymphocyte ratio, PMC4183848, N=49 pts.
        exhausted_total_ratio_mean = 0.366;
        exhausted_total_ratio_std = 0.157;
        % -- measure actual values
        stroma_tumor_ratio_obs =  sum(mySystem.grid.Lf(:)) / numel(mySystem.TU.TUcells);
        lympho_tumor_ratio_obs =  numel(mySystem.IM.IMcells) / numel(mySystem.TU.TUcells);
        macro_tumor_ratio_obs =  numel(mySystem.MP.MPcells) / numel(mySystem.TU.TUcells);
        exhausted_total_ratio_obs = sum(mySystem.IM.IMprop.Kcap == 0) / numel(mySystem.IM.IMcells);
        % -- check deviation from desired values
        p_loss(1) = (stroma_tumor_ratio_obs-stroma_tumor_ratio_mean) / stroma_tumor_ratio_std;
        p_loss(2) = (lympho_tumor_ratio_obs-lympho_tumor_ratio_mean) / lympho_tumor_ratio_std;   
        p_loss(3) = (macro_tumor_ratio_obs-macro_tumor_ratio_mean) / macro_tumor_ratio_std;
        p_loss(4) = (exhausted_total_ratio_obs-exhausted_total_ratio_mean) / exhausted_total_ratio_std;
        
        case 'default_stem' % this includes tumor stem cell fraction, should be used for 3D
        
        % -- define desired values
        if cnst.verbose
            disp('calculating default_stem loss function');
        end
        
        % stroma-tumor ratio, from own data, N=14 patients, 16 Feb 2018
        stroma_tumor_ratio_mean = 0.2263;
        stroma_tumor_ratio_std  = 0.1562;
        % lymphocyte-tumor ratio, from own data, N=14 patients, 16 Feb 2018
        lympho_tumor_ratio_mean = 0.0332;
        lympho_tumor_ratio_std  = 0.0411;
        % macrophage-tumor ratio, from own data, N=14 patients, 16 Feb 2018
        macro_tumor_ratio_mean  = 0.0978;
        macro_tumor_ratio_std   = 0.0825;

        % exhausted-total lymphocyte ratio, PMC4183848, N=49 pts.
        exhausted_total_ratio_mean = 0.366;
        exhausted_total_ratio_std = 0.157;
        % tumor stem cell fraction, 10.1038/nature21713, Extended Data Figure 2
        tumor_stem_cell_fraction_mean = 0.250;
        tumor_stem_cell_fraction_std = 0.200; % not given in the paper -> guessed

        % -- measure actual values
        stroma_tumor_ratio_obs =  sum(mySystem.grid.Lf(:)) / numel(mySystem.TU.TUcells);
        lympho_tumor_ratio_obs =  numel(mySystem.IM.IMcells) / numel(mySystem.TU.TUcells);
        macro_tumor_ratio_obs =  numel(mySystem.MP.MPcells) / numel(mySystem.TU.TUcells);
        exhausted_total_ratio_obs = sum(mySystem.IM.IMprop.Kcap == 0) / numel(mySystem.IM.IMcells);
        tumor_stem_cell_fraction_obs = sum(mySystem.TU.TUprop.isStem == 1) / numel(mySystem.TU.TUcells);

        % -- check deviation from desired values
        p_loss(1) = (stroma_tumor_ratio_obs-stroma_tumor_ratio_mean) / stroma_tumor_ratio_std;
        p_loss(2) = (lympho_tumor_ratio_obs-lympho_tumor_ratio_mean) / lympho_tumor_ratio_std;   
        p_loss(3) = (macro_tumor_ratio_obs-macro_tumor_ratio_mean) / macro_tumor_ratio_std;
        p_loss(4) = (exhausted_total_ratio_obs-exhausted_total_ratio_mean) / exhausted_total_ratio_std;
        p_loss(5) = (tumor_stem_cell_fraction_obs-tumor_stem_cell_fraction_mean) / tumor_stem_cell_fraction_std;

        case 'patient_HEIDEL'
            % load individual patient data and use it for loss function
            
            % stroma-tumor ratio measured in the patient
            stroma_tumor_ratio_mean = cnst.patient.stroma_tumor_ratio;
            stroma_tumor_ratio_std  = cnst.patient.stroma_tumor_ratio * 0.25;
            % lymphocyte-tumor ratio measured in the patient
            lympho_tumor_ratio_mean = cnst.patient.lympho_tumor_ratio;
            lympho_tumor_ratio_std  = cnst.patient.lympho_tumor_ratio * 0.25;
            % macrophage-tumor ratio measured in the patient
            macro_tumor_ratio_mean  = cnst.patient.macro_tumor_ratio;
            macro_tumor_ratio_std   = cnst.patient.macro_tumor_ratio * 0.25;
                      
             % exhausted-total lymphocyte ratio, PMC4183848, N=49 pts.
             exhausted_total_ratio_mean = 0.366;
             exhausted_total_ratio_std = 0.157;
             % tumor stem cell fraction, 10.1038/nature21713, Extended Data Figure 2
             tumor_stem_cell_fraction_mean = 0.250;
             tumor_stem_cell_fraction_std = 0.200;

            % -- measure actual values
            stroma_tumor_ratio_obs =  sum(mySystem.grid.Lf(:)) / numel(mySystem.TU.TUcells);
            lympho_tumor_ratio_obs =  numel(mySystem.IM.IMcells) / numel(mySystem.TU.TUcells);
            macro_tumor_ratio_obs =  numel(mySystem.MP.MPcells) / numel(mySystem.TU.TUcells);
             exhausted_total_ratio_obs = sum(mySystem.IM.IMprop.Kcap == 0) / numel(mySystem.IM.IMcells);
             tumor_stem_cell_fraction_obs = sum(mySystem.TU.TUprop.isStem == 1) / numel(mySystem.TU.TUcells);

            % -- check deviation from desired values
            p_loss(1) = (stroma_tumor_ratio_obs-stroma_tumor_ratio_mean) / stroma_tumor_ratio_std;
            p_loss(2) = (lympho_tumor_ratio_obs-lympho_tumor_ratio_mean) / lympho_tumor_ratio_std;   
            p_loss(3) = (macro_tumor_ratio_obs-macro_tumor_ratio_mean) / macro_tumor_ratio_std;
             p_loss(4) = (exhausted_total_ratio_obs-exhausted_total_ratio_mean) / exhausted_total_ratio_std;
             p_loss(5) = (tumor_stem_cell_fraction_obs-tumor_stem_cell_fraction_mean) / tumor_stem_cell_fraction_std;
        
        case 'antiadju'
             
            % stroma-tumor ratio measured in the patient
            targetLowAnti_mean = cnst.patient.lowAnti;
            targetLowAnti_std  = 0.05;
            % lymphocyte-tumor ratio measured in the patient
            targetLowAdju_mean = cnst.patient.lowAdju;
            targetLowAdju_std  = 0.05;

            % -- measure actual values
            adjuvanticities =  mySystem.grid.AdjuMap(mySystem.TU.TUcells); % get all tumor cell adjuvanticities
            LowAnti_observed  = sum(mySystem.TU.TUprop.Antigen <= mySystem.params.antiThresh) / numel(mySystem.TU.TUcells);
            LowAdju_observed = sum( adjuvanticities <= mySystem.params.adjuThresh)  / numel(mySystem.TU.TUcells);

            % -- check deviation from desired values
            p_loss(1) = (LowAnti_observed-targetLowAnti_mean) / targetLowAnti_std;
            p_loss(2) = (LowAdju_observed-targetLowAdju_mean) / targetLowAdju_std;              
            
        case 'none' % if no loss requested, then no loss returned
            p_loss = []; 
        otherwise
            warning('invalid loss function specified');
    end
    
    % common post-processing
    if cnst.averageOut
     m_loss = sqrt(mean((p_loss).^2)); %total loss is ROOT MEAN SQUARE ERROR
    else 
     m_loss = p_loss;
    end
   %  m_loss = mean(abs(p_loss));
     
    % penalize tumor eradication and NaN or Inf loss 
    if numel(mySystem.TU.TUcells) == 0 || isnan(mean(m_loss)) || isinf(mean(m_loss))
        m_loss = cnst.penalty;
    end
   % mean of partial losses
   if cnst.verbose
    disp(['partial losses: ', num2str(p_loss), ', total loss: ', num2str(m_loss)]);
   end 
else % if no loss requested, then no loss returned
    p_loss = []; m_loss = [];
end

end