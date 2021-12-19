function contentInPatch = PatchContent(ImagesIn,imageColOrder,LEDsUsed,ROIList,mask)
nImgs = size(ImagesIn,3);
if nargin==5 && ~isempty(mask)
    contentMask = mask;
else
    LEDsUsedIdx = [];
    for nn = 1:nImgs
        [ky,kx] = find(imageColOrder == nn);
        if LEDsUsed(ky,kx) == 1
            LEDsUsedIdx = [LEDsUsedIdx,nn];
        end
    end
    onAxisIncoherent = sum(ImagesIn(:,:,LEDsUsedIdx),3);
    [~,threshold] = edge(onAxisIncoherent,'sobel');
    contentMask = edge(onAxisIncoherent,'sobel',threshold * 0.5);
end
contentInPatch = zeros([size(ROIList,1),1]);
for i=1:size(ROIList,1)
    coords = ROIList(i,:);
    if nnz(contentMask(coords(1):coords(2),coords(3):coords(4)))
        contentInPatch(i) = 1;
    end
end