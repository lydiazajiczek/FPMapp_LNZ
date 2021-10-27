function [ImagesIn] = BackgroundRemovingROI(ImagesIn, bck, scaleFactor, LEDs, imageColOrder, systemSetup)
% Function that removes background from darkfield images using 
% mean value of a user selected ROI
%   Inputs:
%       ImagesIn - vector of input images
%       ROI_bg - ROI of background area selected by user
%       LEDS - LEDs used during acquisition
%       imageColOrder - collection order of images
%       systeSetup - system parameters
%   Outputs:
%       I - vector of background removed input images (darkfield only)
%       bck - vector of background values (zero for brightfield images)

ImagesIn = double(ImagesIn);
bck = double(bck);
nImgs = size(ImagesIn,3);
[dF,~,~,~] = Used_LEDs(LEDs,systemSetup,3);
dF = logical(dF);
%dfImageIndices = imageColOrder.*dF;
%idx = dfImageIndices > 0;
%dfImageList(:) = sort(dfImageIndices(idx))

if nImgs == max(max(imageColOrder))
    for i=1:nImgs
        [j,k] = find(imageColOrder==i);
        if(dF(j,k))
            I = ImagesIn(:,:,i);
            bG = bck(i);
            ImagesIn(:,:,i) = I - scaleFactor*bG;
        end
    end
else
    error('Size of input image stack not equal to collection order list')
end

