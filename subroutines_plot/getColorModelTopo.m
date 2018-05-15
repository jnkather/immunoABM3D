
function myColVector = getColorModelTopo(currData,colornorm)
   
    if isnan(currData) 
        myColVector = [0 0 0]; 
    else
    if currData>colornorm, currData =0; end
     
    nlevels = 255;
    currColor = round(currData/colornorm*(nlevels-1))+1;
        
    allColors = redblu(nlevels);
    myColVector = allColors(currColor,:);
    
    end
end

