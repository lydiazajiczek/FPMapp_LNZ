function incoh_img = IncoherentImage(ImagesIn,LEDsUsed,imageColOrder)
%INCOHERENTIMAGE a function to produce an incoherent (sum) image of all
% images chosen for reconstruction.
LEDIdx = sort(nonzeros(imageColOrder.*LEDsUsed));
max_val = double(max(max(max(ImagesIn(:,:,LEDIdx)))));
incoh_img = double(sum(uint32(ImagesIn(:,:,LEDIdx)),3));
incoh_img = incoh_img.*(max_val/max(max(incoh_img)));
end