%%%%function to run FPM reconstruction using set parameters in 
%%%%initialization and initialization2.mat
%%%%required inputs: 
%%%%-filepath to either directory containing TIFFs or single multipage TIFF
%%%%-type of reconstruction ('single_directory', 'single_TIFF', 
%%%%'multi_TIFFs')
%%%%-type of stitching ('single_ROI', 'test_ROIs', 'all_ROIs')
%%%%optional inputs: 
%%%%-'ROILength' (for single_ROI, square works best so only one side, default 256)
%%%%-'ROIVertex' (for single_ROI, otherwise finds central ROI: [x0, y0])
%%%%-'overlap' (for all_ROIs, default is 10%, provide in decimal)
%%%%-'CorrLEDPosns' (filepath to corrected LED positions .mat file, 
%%%%default hardcoded in file)
%%%%-'Keyword' to filter multipage TIFFs to process in folder e.g. 'blue'

function FPM_run(filepath,reconType,stitchType,varargin)

expectedReconTypes = {'single_directory', 'single_TIFF', 'multi_TIFF'};
expectedStitchTypes = {'single_ROI', 'test_ROI', 'all_ROIs'};
defaultOverlap = 0.1;
defaultROILength = 256;
defaultROIVertex = [-1, -1]; %find central ROI if not specified
%defaultCorrLEDPosns = 'C:\Datasets\FPM\11\507983\correction_results.mat';
defaultCorrLEDPosns = 'F:\FPM images\2021\11\508993\leds.mat';

p = inputParser;
validFilepath = @(x) exist(x,'file')>0;
validReconTypes = @(x) any(validatestring(x,expectedReconTypes));
validStitchTypes = @(x) any(validatestring(x,expectedStitchTypes));
validOverlap = @(x) isnumeric(x) && (x >=0 ) && (x < 1);
validROILength = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validROIVertex = @(x) isvector(x) && length(x) == 2 ...
                        && all(x > 0) && all(rem(x,1)==0);
validKeyword = @(x) isstring(x) || ischar(x);
addRequired(p,'filepath',validFilepath)
addRequired(p,'reconType',validReconTypes)
addRequired(p,'stitchType',validStitchTypes)
addParameter(p,'overlap',defaultOverlap,validOverlap);
addParameter(p,'ROILength',defaultROILength,validROILength);
addParameter(p,'ROIVertex',defaultROIVertex,validROIVertex);
addParameter(p,'CorrLEDPosns',defaultCorrLEDPosns,validFilepath)
addParameter(p,'Keyword','',validKeyword);
parse(p,filepath,reconType,stitchType,varargin{:});

overlap = p.Results.overlap;
sx = p.Results.ROILength;
ROIVertex = p.Results.ROIVertex;
keyword = p.Results.Keyword;

%%%%Initialization
%init_path = 'C:\Users\lydia\Code\github_repos\FPMapp_LNZ\';
init_path = 'F:\FPMapp_LNZ\';
addpath([init_path '\Algorithm functions']);
addpath([init_path '\Algorithm functions\BM3D']);
addpath([init_path '\Algorithm functions\sort_nat']);
addpath([init_path '\Algorithm functions\DPC']);
addpath([init_path '\Algorithm functions\DPC\dpc_functions']);
addpath([init_path '\Algorithm functions\minFunc']);
addpath([init_path '\GUI functions']);

load([init_path '/settings/settings.mat'], ...
    'imageColOrder', 'LEDs', 'LEDsUsed', 'options', 'systemSetup'); 
load([init_path '/settings/settings2.mat'], 'others'); 

options.dfBackground = false; %can change this obviously
options.LEDcorrection = '0';
options.maxIter = 5;
tmp0 = load(p.Results.CorrLEDPosns,'svdIdx');
others.showIterResult = false;
others.loadPrevLEDPos = true;
others.saveIterations = false;
%others.iterationsSaveDir = 'C:\Datasets\FPM\10\08\513282\green_1_reconstruction\';

FOVList = [];
switch reconType
    case 'single_directory'
        if ~isfolder(filepath)
            error('Filepath is not a directory.')
        end
        FOVList = [FOVList string(filepath)]; %#ok<*NBRAK>
%         saveDir = filepath;%[filepath '\reconstruction\'];
    case 'single_TIFF'
        if isfolder(filepath)
            error('Filepath is a directory, not a file.')
        end
        FOVList = [FOVList string(filepath)];
        [folder, ~, ~] = fileparts(filepath);
%         saveDir = folder;
%         saveDir = [folder '\' filename '_reconstruction\'];
    case 'multi_TIFF'
        if ~isfolder(filepath)
            error ('Filepath is not a directory.')
        end
        FOVList = LoadDir(filepath,keyword);
%         saveDir = filepath;
        if isempty(FOVList)
            error('Directory contains no multipage TIFFs.')
        end
end

%%%%Loop through all FOVs
for i=1:length(FOVList)
    switch reconType
        case 'single_directory'
            folder = string(FOVList(i));
            [I,~] = LoadImages(folder);
            filename = '';
        case 'single_TIFF'
            I = LoadImage(string(FOVList(i)));
            [folder, filename, ~] = fileparts(FOVList(i));
        case 'multi_TIFF'
            
            [folder, filename, ~] = fileparts(FOVList(i));
            split_filename = split(filename,'_');
            colour = split_filename(1);
            img_num = split_filename(2);
        
            img_num = str2num(img_num);
            if img_num < 14 || img_num > 126
                continue;
            end
        
            switch lower(colour)
                case 'red'
                    systemSetup.lambda = 0.6292;
                case 'green'
                    systemSetup.lambda = 0.53;
                case 'blue'
                    systemSetup.lambda = 0.475;
            end
            I = LoadImage(FOVList(i));
%             saveDir = [folder '\' filename '_reconstruction\'];
    end

    
    
    [ny,nx,~] = size(I);
    
    %%%%obtain summed incoherent image
    incoh_img = sum(I,3);
    incoh_img = uint16(incoh_img.*double(max(max(max(I))))/max(max(incoh_img)));
     
    %%%%find all ROIs including max contrast ROI
    contrast_max = 0;
    ROI_max = [];
    j = 1;
    ix = 1;
    iy = 1;
    while (nx-ix+1) > sx  
        while (ny-iy+1) > sx
            ROI = [iy,iy+sx-1,ix,ix+sx-1];
            patch = imcrop(incoh_img,ROI);
            max_patch = max(max(patch));
            min_patch = min(min(patch));
            %contrast = double(max_patch-min_patch)/double(max_patch+min_patch);
            contrast = std2(patch);
            if contrast > contrast_max
                contrast_max = contrast;
                ROI_max = ROI;
            end
            ROIList(j,:) = [iy,iy+sx-1,ix,ix+sx-1]; %#ok<AGROW>
            iy = (iy+sx-overlap*sx);  
        j = j+1;
        end
        ix = (ix+sx-overlap*sx); 
        iy = 1;
    end
    
    if ~isempty(filename)
        filename = strcat(filename,'_');
    end
        
    switch stitchType
        case 'single_ROI'
            if any(ROIVertex < 0)
                %%%%find central ROI
                ROI = centralROI(nx,ny,sx,sx);
                disp(['reconstructing central ROI: [' num2str(ROI(1)) ' ' ...
                     num2str(ROI(2)) ' ' num2str(ROI(3)) ' ' num2str(ROI(4)) ']'])
            else
                ROI(1) = ROIVertex(1);
                ROI(2) = ROIVertex(2);
                ROI(3) = sx;
                ROI(4) = sx;
            end
            %ROIList = [ROI];
            [rec_object,phase,rec_pupil,~,~,~,~,~] = AlgorithmManual(...
                I, LEDs, imageColOrder, LEDsUsed, ROI, [], ...
                systemSetup, options, others, tmp0.svdIdx, 2);
            imwrite_uint16(imcrop(incoh_img,ROI),strcat(folder,'\',filename,'inc_sum.tif'))
            imwrite_uint16(abs(rec_object),strcat(folder,'\',filename,'amplitude.tif'));
            imwrite_float(single(phase),strcat(folder,'\',filename,'phase.tif'));
            imwrite_float(single(angle(rec_pupil)),strcat(folder,'\',filename,'pupil_phase.tif'));
            
        case 'test_ROI'
%             ROIList = [ROI_max];
            ROI(1) = uint16(ROI_max(3));
            ROI(2) = uint16(ROI_max(1));
            ROI(3) = sx;
            ROI(4) = sx;
            
            writematrix(ROI,strcat(folder,'\',filename,'ROI.txt'),'Delimiter',',');
            
            disp(['reconstructing max contrast ROI: [' num2str(ROI(1)) ' ' ...
                 num2str(ROI(2)) ' ' num2str(ROI(3)) ' ' num2str(ROI(4)) ']'])
            [rec_object,phase,rec_pupil,~,~,~,~,~] = AlgorithmManual(...
                I, LEDs, imageColOrder, LEDsUsed, ROI, [], ...
                systemSetup, options, others, tmp0.svdIdx, 2);
            
            imwrite_uint16(imcrop(incoh_img,ROI),strcat(folder,'\',filename,'test_inc_sum.tif'))
            imwrite_uint16(abs(rec_object),strcat(folder,'\',filename,'test_amplitude.tif'));
            imwrite_float(single(phase),strcat(folder,'\',filename,'test_phase.tif'));
            imwrite_float(single(angle(rec_pupil)),strcat(folder,'\',filename,'test_pupil_phase.tif'));
            
            %I, LEDs, imageColOrder, LEDsUsed, ROIList(90,:), ...
        case 'all_ROIs'
            saveDir = strcat(folder,'\',filename,'reconstruction\');
            %%%%make save directory
            if exist(saveDir,'dir')==0
                mkdir(saveDir)
            end
            AlgorithmManualStitch(...
                        I, LEDs, imageColOrder, LEDsUsed, ROIList, ...
                        [], saveDir, systemSetup, options, ...
                        others, tmp0.svdIdx, 2);
            
    end
            %ROIList = zeros(5,4);
            %%%%top left
            %ROIList(1,:) = [1, 1, sx, sx];
            %%%%top right
            %ROIList(2,:) = [1, nx-sx+1, sx, sx];
            %%%%central ROI
            %ROIList(3,:) = centralROI(nx,ny,sx,sx);
            %%%%bottom left
            %ROIList(4,:) = [ny-sx+1, 1, sx, sx]; 
            %%%%bottom right
            %ROIList(5,:) = [ny-sx+1, nx-sx+1, sx, sx];   
end
end

function ROI = centralROI(nx,ny,sx,sy)
    ROI = zeros(4,1);
    ROI(1) = floor(ny/2) - ceil(sy/2);
    ROI(2) = floor(nx/2) - ceil(sx/2);
    ROI(3) = sy;
    ROI(4) = sx;
end
