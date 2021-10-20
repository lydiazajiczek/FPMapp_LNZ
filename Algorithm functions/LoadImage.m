function I = LoadImage(filePath)
% Function that loads large TIFF file
%   Inputs:
%       filePath - path of TIFF file
%   Outputs:
%       I - loaded images (3D matrix)   

f = waitbar(0,'Loading big TIFF ... ');
I = double(TIFFStack(filePath));
try
    close(f)
end
