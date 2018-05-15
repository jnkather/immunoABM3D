% JN Kather Feb 2018

function [sysOut, summaryPreTherapy, summaryEnd] = ...
    runWithTreatment(sysTempl, cnst, therapy, masterExpname, numexp)

saveImage = cnst.saveImage;
saveVideo = cnst.saveVideo;

% prepare system package for each run
for i=1:numexp
    sysTempl.params.initialSeed = NaN;
    allSys{i}= sysTempl;
end

%% run experiments
for i=1:numexp
    if cnst.verbose
    disp(['starting experimental run ',num2str(i)]);
    end
    
    currSys = allSys{i}; % retrieve system template
    expname = [masterExpname,'_iter',num2str(i)];    % experiment name that will be used to save results

    if therapy.active % apply therapy
    sysOut{i} = currSys;
    cnst_override = cnst;
        for j = 1:numel(therapy.time) % for each therapy breakpoint
            if j == 1 % first therapy breakpoint
                % run the system - PRE THERAPY
                disp('starting PRE-therapy run...');
                cnst_override.nSteps = therapy.time(j); % run up to the first therapy time point
                disp(['will run for ',num2str(cnst_override.nSteps),' iterations']);
                cnst_override.createNewSystem = true;
                [sysOut{i}, ~, summaryPreTherapy{i}, imWin, ~, ~] = ...
                 runSystem(sysOut{i},cnst_override,expname,saveImage,saveVideo); 
                if ~imWin % if tumor still alive
                disp('starting POST-first-therapy run...');
                %disp('press any key to continue...');pause
                % apply first therapy
                sysOut{i} = executeTherapy(sysOut{i},therapy.drug{j},therapy.strength(j),cnst.verbose);
                % run the system - POST THERAPY
                if j+1 <= numel(therapy.time) % if there will be another therapy
                    cnst_override.nSteps = therapy.time(j+1) - therapy.time(j);
                else
                    cnst_override.nSteps = cnst.nSteps - therapy.time(j);
                end
                disp(['will run for ',num2str(cnst_override.nSteps),' iterations']);
                cnst_override.requireAlive = 1; % added march 2018
                cnst_override.createNewSystem = false;
                cnst_override.iterationsAdd = therapy.time(j);
                [sysOut{i}, ~, summaryEnd{i}, ~, ~, ~] = ...
                 runSystem(sysOut{i},cnst_override,expname,saveImage,saveVideo);  
                else
                    error('tumor was eradicated');
                end
            elseif j==numel(therapy.time) % if this is the last therapy
                % apply last therapy
                disp('starting the POST-LAST-therapy run...');
                sysOut{i} = executeTherapy(sysOut{i},therapy.drug{j},therapy.strength(j));
                % run the system - POST THERAPY
                cnst_override.nSteps = cnst.nSteps - therapy.time(j);
                disp(['will run for ',num2str(cnst_override.nSteps),' iterations']);
                cnst_override.requireAlive = 1; % added march 2018
                cnst_override.createNewSystem = false;
                cnst_override.iterationsAdd = therapy.time(j);
                [sysOut{i}, ~, summaryEnd{i}, ~, ~, ~] = ...
                 runSystem(sysOut{i},cnst_override,expname,saveImage,saveVideo);    
            else % if this is not the first and not the last
                % apply last therapy
                disp('starting the POST-jth-therapy run...');
                sysOut{i} = executeTherapy(sysOut{i},therapy.drug{j},therapy.strength(j),cnst.verbose);
                % run the system - POST THERAPY
                cnst_override.nSteps = therapy.time(j+1) - therapy.time(j);
                disp(['will run for ',num2str(cnst_override.nSteps),' iterations']);
                cnst_override.requireAlive = 1; % added march 2018
                cnst_override.createNewSystem = false;
                cnst_override.iterationsAdd = therapy.time(j);
                [sysOut{i}, ~, ~, ~, ~, ~] = ...
                 runSystem(sysOut{i},cnst_override,expname,saveImage,saveVideo);    
            end
        end
    else % do not apply therapy
    [sysOut{i}, ~, summaryEnd{i}, ~, ~, ~] = ...
        runSystem(currSys,cnst,expname,saveImage,saveVideo); % run the system
    summaryPreTherapy = summaryEnd; % if no therapy, no difference here
    end

 %   paramsOut{i} = sysOut{i}.params;
end


end