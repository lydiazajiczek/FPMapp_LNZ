function contentInPatch = PatchContent(ImagesIn,imageColOrder,LEDsUsed,ROIList)
nImgs = size(ImagesIn,3);
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
se90 = strel('line',3,90);
se0 = strel('line',3,0);
%contentMask = imdilate(contentMask,[se90 se0]);
contentInPatch = zeros([nImgs,1]);
for i=1:size(ROIList,1)
    coords = ROIList(i,:);
    if nnz(contentMask(coords(1):coords(2),coords(3):coords(4)))
        contentInPatch(i) = 1;
    end
end