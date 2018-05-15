% Copyright (c) 2017-2018, Jakob Nikolas Kather. 
% 
% Please cite our publication:
% "Large-scale database mining reveals hidden trends and future directions
% for cancer immunotherapy", DOI 10.1080/2162402X.2018.1444412
% 
% License: please refer to the license file in the root directory
%
% -------------------------------------------------------------
%
% this function returns up to eight colors from the color brewer color
% maps, plese refer to http://colorbrewer2.org/
% 

function cmapout = brewer1(numC)

hexout = {'#fb6a4a', '#de2d26', '#a50f15', ...
'#6baed6','#3182bd','#08519c','#54278f',...
'#31a354','#fe9929','#000000'};

for i=1:numC
   
    currColInd = mod(i-1,numel(hexout))+1;
    cmapout(i,:) = hex2rgb(hexout{currColInd});

end

cmapout = cmapout(1:numC,:);

end