function ROI = ROISelection(in,bg)
% Function that selects Region Of Interest
%   Input:
%       in - image
%       bg - is it background ROI selection?
%   Output:
%       ROI - Region Of Interest       
%           ROI = [x0,y0,xSize,ySize];
figure(100);

if ~bg
    set(gcf,'Name','ROI selection');
else
    set(gcf,'Name','Background ROI selection');
end
set(gcf,'NumberTitle','off');
img = imagesc(in); colormap gray; 
if ~bg
    title('ROI selection');
else
    title('Background ROI selection');
end
[~,~,~,ROI] = imcrop(img);
ROI = round(ROI);
if ~isempty(ROI)
    if ( mod(ROI(3),2) == 1 ), ROI(3) = ROI(3)+1; end 
    if ( mod(ROI(4),2) == 1 ), ROI(4) = ROI(4)+1; end
    delete(gcf)
end

end
