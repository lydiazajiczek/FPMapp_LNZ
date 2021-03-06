%create top, bottom, right images for quantitative DPC
%todo rewrite using indices rather than for loops
function [IDPC, incoh_img] = CreateDPCImgs(ImagesIn,LEDsUsed,imageColOrder,aberrationCorrection)

%%calculate incoherent image here for use as amplitude initialiation
%%and for normalization during stitching

incoh_img = IncoherentImage(ImagesIn,LEDsUsed,imageColOrder);

[nx, ny] = size(LEDsUsed);
IDPC = double(zeros([size(ImagesIn,1) size(ImagesIn,2) 4]));

%top image
idx = LEDsUsed(1:nx/2,:) > 0;
colOrder = imageColOrder(1:nx/2,:);
index = colOrder(idx);
maxVal = 0;
for i=1:length(index)
    IDPC(:,:,1) = IDPC(:,:,1) + double(ImagesIn(:,:,index(i)));
    if max(max(double(ImagesIn(:,:,index(i))))) > maxVal
        maxVal = max(max(double(ImagesIn(:,:,index(i)))));
    end
end
IDPC(:,:,1) = IDPC(:,:,1).*maxVal/(max(max(IDPC(:,:,1))));
clear index

%bottom image
idx = LEDsUsed((nx/2+1):end,:) > 0;
colOrder = imageColOrder((nx/2+1):end,:);
index = colOrder(idx);
maxVal = 0;
for i=1:length(index)
    IDPC(:,:,2) = IDPC(:,:,2) + double(ImagesIn(:,:,index(i)));
    if max(max(double(ImagesIn(:,:,index(i))))) > maxVal
        maxVal = max(max(double(ImagesIn(:,:,index(i)))));
    end
end
IDPC(:,:,2) = IDPC(:,:,2).*maxVal/(max(max(IDPC(:,:,2))));
clear index

%right image
idx = LEDsUsed(:,(nx/2+1):end) > 0;
colOrder = imageColOrder(:,(ny/2+1):end);
index = colOrder(idx);
maxVal = 0;
for i=1:length(index)
    IDPC(:,:,3) = IDPC(:,:,3) + double(ImagesIn(:,:,index(i)));
    if max(max(double(ImagesIn(:,:,index(i))))) > maxVal
        maxVal = max(max(double(ImagesIn(:,:,index(i)))));
    end
end
IDPC(:,:,3) = IDPC(:,:,3).*maxVal/(max(max(IDPC(:,:,3))));
clear index

if ~aberrationCorrection
    %left image
    idx = LEDsUsed(:,1:nx/2) > 0;
    colOrder = imageColOrder(:,1:ny/2);
    index = colOrder(idx);
    maxVal = 0;
    for i=1:length(index)
        IDPC(:,:,4) = IDPC(:,:,4) + double(ImagesIn(:,:,index(i)));
        if max(max(double(ImagesIn(:,:,index(i))))) > maxVal
            maxVal = max(max(double(ImagesIn(:,:,index(i)))));
        end
    end
    IDPC(:,:,4) = IDPC(:,:,4).*maxVal/(max(max(IDPC(:,:,4))));
else
    %central coherent image
    IDPC(:,:,4) = ImagesIn(:,:,ceil(size(ImagesIn,3)/2));
end
