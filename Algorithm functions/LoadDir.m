function FOVList = LoadDir(loadDirectory)
% Function that checks a directory for big TIFFs
%   Inputs:
%       loadDirectory - directory where are input FPM images
%   Outputs:
%       FOVList - list of paths to big TIFFs to be loaded with LoadImage.m

fileList = dir(loadDirectory);
fileList(2) = []; fileList(1) = [];
a = string(zeros(length(fileList),1));
for k = 1:length(fileList)
    a(k) = fileList(k).name;
end
fileList = sort_nat(a);
clear k a;

%%%filter file list to just contain TIFF files
idx = contains(fileList,'.tif');
imageList = fileList(idx);
clear fileList;

%%%now check to make sure they are big TIFFs
FOVList = [];
for i=1:length(imageList)
    filename_full = strcat(loadDirectory,'\',imageList(i));
    if length(imfinfo(filename_full)) > 1
        FOVList = [FOVList filename_full]; %#ok<AGROW>
    end
end