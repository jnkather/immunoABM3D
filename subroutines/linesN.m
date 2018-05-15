% color map

function cmap = linesN(fh,numC,reps)

     cmap= fh(numC); 
     cmap(end,:) = 0; % black as last color
     replines = ones(numC,1)*reps;
     indices = accumarray(cumsum([1;replines]),1);
     cmap = cmap(cumsum(indices(1:end-1)),:);
    
     
end