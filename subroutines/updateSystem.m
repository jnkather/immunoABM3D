
function [mySystem,currentSummary] = updateSystem(mySystem,state,i,cnst)

    % --- copy all variables back to mySystem ---
    % tumor cell coordinates and properties ------------------------------
    mySystem.TU.TUcells = state.TUcells.position;         % list of tumor cells
    mySystem.TU.TUprop.isStem = state.TUcells.isStem;      % tumor cell stemness
    mySystem.TU.TUprop.Pcap = state.TUcells.Pcap;          % remaining cell divs
    mySystem.TU.TUprop.Antigen = state.TUcells.Antigen;    % cell antigenicity
    mySystem.TU.TUprop.damage = state.TUcells.damage;    % cell antigenicity
    % lymphocyte coordinates and properties ------------------------------
    mySystem.IM.IMcells = state.Lymphocytes.position;                  % list of lymphocytes
    mySystem.IM.IMprop.Kcap = state.Lymphocytes.Kcap;          % remaining kills
    mySystem.IM.IMprop.Pcap = state.Lymphocytes.Pcap;          % remaining cell divs
    mySystem.IM.IMprop.engaged = state.Lymphocytes.engaged;    % lymphocyte engagement
    % macrophage coordinates and properties ------------------------------
    mySystem.MP.MPcells = state.Macrophages.position;                  % list of macrophages
    mySystem.MP.MPprop.Pcap = state.Macrophages.Pcap;          % proliferation capac.
    mySystem.MP.MPprop.State = state.Macrophages.state;        % differentiation state
    % global system maps -------------------------------------------------
    mySystem.grid.ChtaxMap = state.env.ChtaxMap;              % chemotaxis map
    mySystem.grid.HypoxMap = state.env.NecroMap;              % hypoxia map
    mySystem.grid.AdjuMap = state.env.AdjuMap;                % adjuvanticity map
    % necrosis and fibrosis ----------------------------------------------
    mySystem.grid.Ln = state.env.Ln;                          % necrosis map
    mySystem.grid.Lf = state.env.Lf;                          % fibrosis (stroma) map
    % internal grids for faster computation of occupancy -----------------
    mySystem.grid.L = state.env.L;                            % global occupancy map
    mySystem.grid.Lt = false(size(state.env.L));
    mySystem.grid.Lt(mySystem.TU.TUcells) = true;                          % tumor cell map
    mySystem.grid.StepsDone = i;    

    % create immune grid
    mySystem.grid.Li = false(size(state.env.L));
    mySystem.grid.Li(mySystem.IM.IMcells) = true;
    
    if (mySystem.is3D) % summarize system - 3D
        currentSummary = summarizeSystem_3D(mySystem,cnst);
    else               % summarize system - 2D
        currentSummary = summarizeSystem_2D(mySystem,cnst);
    end
    
end
