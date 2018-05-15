function sysOut = executeTherapy(sysOut,myDrug,strength,verbose)

% strength varies between 0 and 20
if verbose
disp(['executing therapy, strength is ', num2str(strength)]);
end

for i = 1:numel(myDrug)
    disp(['     adding drug #',num2str(i)]);
    switch myDrug{i}
        case 'CTRL'
            disp('     --- no change ---');
        case 'CH'
            disp('     adding chemotherapy (increase cell death probability upon proliferation attempt)'); 
            % arbitrary constant factor of 0.01
            sysOut.params.CellDieAtProl = 0.012 * strength; 
        case 'PD'
            disp('     adding anti-pd1 (increase exhaustion limit + un-exhaust cells)');  % Ref. Chen&Mellman
            sysOut.params.IMkmax = min(round(sysOut.params.IMkmax * 2 * strength),250); 
            sysOut.IM.IMprop.Kcap = sysOut.params.IMkmax*ones(size(sysOut.IM.IMprop.Kcap)); 
            
        case 'RE'  % repolarization agent, ref. Halama et al. Cancer Cell 2016
            % constant factor of 3.3 so that at strength 20, MPrepola is 1
            disp('     adding repolarization agent (increase macrophage repolarization)'); 
            sysOut.params.MPprepola = sysOut.params.MPprepola * (6/20) * strength;
        case 'SP'
            disp('     adding permeabilization agent (permeabilize stroma)'); 
            sysOut.params.stromaPerm = strength * 0.04; % should be 0.8 at max, linear
        case 'LB' % lymphocyte boost
            disp('     performing lymphocyte boost (increase influx)'); 
            sysOut.params.IMrateDynamic = sysOut.params.IMrateDynamic * strength;      
        case 'MD' % macrophage depletion
            disp('     remove all existing macrophages once'); 
            sysOut.MP.MPcells = []; 
            sysOut.MP.MPprop.Pcap = [];
            sysOut.MP.MPprop.State = [];        
        case 'MI' % movement inhibition
            disp('     decreasing tumor cell movement'); 
            sysOut.params.TUpmig = sysOut.params.TUpmig - sysOut.params.TUpmig * strength/20; 
        otherwise
            error('invalid therapy type');
    end
            
    
end

end