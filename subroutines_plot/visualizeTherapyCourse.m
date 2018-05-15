function visualizeTherapyCourse(masterExpname,therapy,cnst)

yoffset = 0.15;
xoffset = 25;
    figure()
    set(gcf,'Color','w');
    suptitle(strrep(masterExpname,'_',' '));
    title('treatment course');
    hold on
    
    % plot timeline
    plot([0 cnst.nSteps],[0 0],'k-','LineWidth',2);
    % plot drawWhen points
    xdrawWhen = 0:cnst.drawWhen:cnst.nSteps;
    ydrawWhen = xdrawWhen*0;
    plot(xdrawWhen,ydrawWhen,'kx','LineWidth',2,'MarkerSize',10);
    % label drawWhen points
    text(xdrawWhen,ydrawWhen+yoffset,num2str((xdrawWhen/2)'),'HorizontalAlignment','center');
    text(max(xdrawWhen)+xoffset,yoffset,'days','FontWeight','bold');
    text(max(xdrawWhen)+xoffset,0,'timeline','FontWeight','bold');
    
    % add therapy
    if therapy.active
        for i = 1:numel(therapy.time) % plot each therapy
            plot(therapy.time(i),-yoffset/3*2,'ro','MarkerSize',8,'LineWidth',2);
           % addyoffset =0;
           % text(therapy.time(i),-yoffset/3*4-addyoffset,therapy.drug{i},'HorizontalAlignment','center');
            for j = 1:numel(therapy.drug{i}) % label therapy
                if j>1,addyoffset = j*0.04;else, addyoffset=0; end
                text(therapy.time(i),-yoffset/3*4-addyoffset,therapy.drug{i}{j},'HorizontalAlignment','center');
            end
        end
        text(max(xdrawWhen)+xoffset,-yoffset/3*4,'treatment','FontWeight','bold');
    else
        disp('therapy is not active... will not plot anything');
    end
    
    axis([0-2*xoffset,max(xdrawWhen)+2*xoffset,-5*yoffset,5*yoffset]);
    axis off
        
end