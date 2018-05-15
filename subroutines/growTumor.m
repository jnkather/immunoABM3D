% JN Kather & J Poleszczuk 2017, jakob.kather@nct-heidelberg.de

function [mySystem, finalImage, finalSummary, imWin, fcount] = ...
    growTumor(mySystem,cnst)
% growTumor_2D performs the actual agent-based modeling in 2D
%       input are two structures: mySystem, defining the initial state of
%       the system; and cnst, defining some hyperparameters

% START PREPARATIONS -------------------------------------------
SIMengine = SIMengine_interface();
if cnst.createNewSystem % create a new (empty) system
    initialize(SIMengine,mySystem,cnst);
else % use existing system and grow it further
    initializeFromState(SIMengine,mySystem,cnst);
end
% END PREPARATIONS -------------------------------------------

% START INITIALIZE AUX VARIABLES  ----------------------------------------
imWin = 0;              % set immune win flag to 0
fcount = 0;             % frame counter for video export
finalImage = [];        % set empty resulting image
% END INITIALIZE AUX VARIABLES   -----------------------------------------

% START ITERATION
for i = 1:cnst.nSteps % iterate through time steps
    % START TUMOR CELL ROUND ------------------------------------------------
    TU_go_grow_die(SIMengine); %performing tumor cell action, dying cells are stored internally
    modulateAdjuvanticity(SIMengine); %dying tumor cells increase adjuvanticity
    % END TUMOR CELL ROUND ---------------------------------------------------
    
    % START MODIFY PARAMETER MAPS --------------------------------------------
    updateNecroMap(SIMengine);
    updateChemoMap(SIMengine);
    decayAdjuvanicityMap(SIMengine);
    % END MODIFY PARAMETER MAPS
    
    % START LYMPHOCYTE AND MACROPHAGE ROUND ------------------------------------------------
    % factor for immune cell influx depends linearly on number of tumor cells

    mySystem.params.IMinfluxModifier = round(mySystem.params.IMrateDynamic * double(TUcellsNum(SIMengine)));
    mySystem.params.MPinfluxModifier = round(mySystem.params.MPrateDynamic * double(TUcellsNum(SIMengine)));
    
     %my proposal
     %if (TUcellsNum(SIMengine) > someThreshold) %only if enough TUcell
     %  mySystem.params.IMinfluxModifier = round(mySystem.params.IMrateDynamic * double(TUcellsNum(SIMengine))/(anotherParam+double(TUcellsNum(SIMengine)));
     %end
     %same for macrophages (?)

    % randomly trigger lymphocyte influx N times depending on currIMinfluxModifier
    if rand() <= mySystem.params.IMinfluxProb
        IMinflux(SIMengine, mySystem.params.IMinfluxModifier);
    end
    % randomly trigger lymphocyte influx N times depending on currMPinfluxModifier
    if rand() <= mySystem.params.MPinfluxProb 
        MPinflux(SIMengine, mySystem.params.MPinfluxModifier);
    end
    
    lymphocytesAct(SIMengine);
    macrophagesAct(SIMengine);
    % END LYMPHOCYTE AND MACROPHAGE ROUND --------------------------------------------------
    
    % START FIBROSIS ------------------------------------------------------
    seedFibrosis(SIMengine);
    % END FIBROSIS ----------------------------------------------------------
    
    % START DRAWING ---------------------------------------------------------
    tumorIsGone = TUcellsNum(SIMengine) == 0;
    if TUcellsNum(SIMengine) >= cnst.maxCells || (mod(i-1,cnst.drawWhen)==cnst.drawWhen-1) || tumorIsGone % plot status after N epochs
        fcount = fcount+1; % video frame counter increases
        if cnst.verbose
            disp(['finished iteration ',num2str(i)]);
        end
        % getting state of the system and export current state of the system
        [mySystem,finalSummary{fcount}] = updateSystem(mySystem,getState(SIMengine),i,cnst);
        if cnst.verbose % enforce drawing and create image from plot
           if mySystem.is3D
                visualize_balls_3D_blank(mySystem,cnst,finalSummary);
            else
                visualize_balls_2D_blank(mySystem,cnst,finalSummary);
            end
            drawnow
            currFrame = getframe(gcf);              % retrieve frame to save later
            finalImage{fcount} = currFrame.cdata;
        else finalImage = [];
        end
        
    end
    % END DRAWING -----------------------------------------------------------
    
    if TUcellsNum(SIMengine) >= cnst.maxCells
       if cnst.verbose
        disp('I hit the limit for max number of tumor cells. will break the loop.'); 
       end
       break; 
    end
    
    % if there are no tumor cells anymore then the game is over
    if tumorIsGone
        if cnst.verbose
            disp('Immune cells win');
        end
        imWin = i;
        return
    end
end % END ITERATION

%CLEARING SIM ENGINE
clear SIMengine;

end % END FUNCTION
