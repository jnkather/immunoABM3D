function writeMyVideo(allImages,finalPath,repeatFrames)

v = VideoWriter(finalPath,'MPEG-4');
open(v)

disp('size of final image container:');
size(allImages)

%targetVideoWidth = 600;
%targetVideoHeight = size(allImages{1},2)* size(allImages{1},1)/targetVideoWidth;

% write all frames to video
for i=1:numel(allImages) % iterate all images
    %disp(['current original frame size = ',num2str(size(allImages{i},1)),'x',...
    %    num2str(size(allImages{i},2)),'... resizing to target size.']);
    %currFrame = imresize(allImages{i},[targetVideoWidth,targetVideoHeight]);
    for j=1:repeatFrames % repeat the frame N times
    writeVideo(v,allImages{i});
    end
end
    
close(v)

end