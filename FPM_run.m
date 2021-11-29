%%%%function to run FPM reconstruction using set parameters in 
%%%%initialization and initialization2.mat
%%%%required inputs: 
%%%%-filepath to either directory containing TIFFs or single multipage TIFF
%%%%-type of reconstruction ('single_directory', 'single_TIFF', 
%%%%'multi_TIFFs')
%%%%-type of stitching ('single_ROI', 'test_ROIs', 'all_ROIs')
%%%%optional inputs: 
%%%%-'ROILength' (square works best so only one parameter, default 256)
%%%%-'ROIVertex' (otherwise finds central ROI)
%%%%-'overlap' (default is 10%, provide in decimal)
%%%%-'CorrLEDPosns' (filepath to corrected LED positions .mat file, 
%%%%default hardcoded in file)

function FPM_run(filepath,reconType,stitchType,varargin)

expectedReconTypes = {'single_directory', 'single_TIFF', 'multi_TIFFs'};
expectedStitchTypes = {'single_ROI', 'test_ROI', 'all_ROIs'};
defaultOverlap = 0.1;
defaultROILength = 256;
defaultROIVertex = [-1, -1]; %find central ROI if not specified
defaultCorrLEDPosns = 'C:\Datasets\FPM\11\507983\correction_results.mat';

p = inputParser;
validFilepath = @(x) exist(x,'file')>0;
validReconTypes = @(x) any(validatestring(x,expectedReconTypes));
validStitchTypes = @(x) any(validatestring(x,expectedStitchTypes));
validOverlap = @(x) isnumeric(x) && (x >=0 ) && (x < 1);
validROILength = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validROIVertex = @(x) isvector(x) && length(x) == 2 && all(x > 0);
addRequired(p,'filepath',validFilepath)
addRequired(p,'reconType',validReconTypes)
addRequired(p,'stitchType',validStitchTypes)
addParameter(p,'overlap',defaultOverlap,validOverlap);
addParameter(p,'ROILength',defaultROILength,validROILength);
addParameter(p,'ROIVertex',defaultROIVertex,validROIVertex);
addParameter(p,'CorrLEDPosns',defaultCorrLEDPosns,validFilepath)
parse(p,filepath,reconType,stitchType,varargin{:});

overlap = p.Results.overlap;
sx = p.Results.ROILength;
ROIVertex = p.Results.ROIVertex;

%%%%Initialization
init_path = 'C:\Users\lydia\Code\github_repos\FPMapp_LNZ\';
addpath([init_path '\Algorithm functions']);
addpath([init_path '\Algorithm functions\BM3D']);
addpath([init_path '\Algorithm functions\sort_nat']);
addpath([init_path '\Algorithm functions\DPC']);
addpath([init_path '\Algorithm functions\DPC\dpc_functions']);
addpath([init_path '\Algorithm functions\minFunc']);
addpath([init_path '\GUI functions']);

load([init_path '/initialization.mat'], ...
    'imageColOrder', 'LEDs', 'LEDsUsed', 'options', 'systemSetup'); 
load([init_path '/initialization2.mat'], 'others'); 

options.dfBackground = false; %can change this obviously
options.LEDcorrection = '0';
tmp0 = load(p.Results.CorrLEDPosns,'svdIdx');
others.showIterResult = false; %doesn't make sense for parfor
others.loadPrevLEDPos = true;
others.saveIterations = true; %see what's going on
others.iterationsSaveDir = 'C:\Datasets\FPM\10\08\513282\green_1_reconstruction\';

FOVList = [];
switch reconType
    case 'single_directory'
        if ~isfolder(filepath)
            error('Filepath is not a directory.')
        end
        FOVList = [FOVList string(filepath)]; %#ok<*NBRAK>
        saveDir = [filepath '\reconstruction\'];
    case 'single_TIFF'
        if isfolder(filepath)
            error('Filepath is a directory, not a file.')
        end
        FOVList = [FOVList string(filepath)];
        [folder, filename, ~] = fileparts(filepath);
        saveDir = [folder '\' filename '_reconstruction\'];
    case 'multi_TIFF'
        if ~isfolder(filepath)
            error ('Filepath is not a directory.')
        end
        FOVList = LoadDir(filepath);
        if isempty(FOVList)
            error('Directory contains no multipage TIFFs.')
        end
end

%%%%Loop through all FOVs
for i=1:length(FOVList)
    switch reconType
        case 'single_directory'
            [I,~] = LoadImages(string(FOVList(i)));
        case 'single_TIFF'
            I = LoadImage(string(FOVList(i)));
        case 'multi_TIFF'
            I = LoadImage(FOVList(i));
            [folder, filename, ~] = fileparts(FOVList(i));
            saveDir = [folder '\' filename '_reconstruction\'];
    end
    %%%%make save directory
    if exist(saveDir,'dir')==0
        mkdir(saveDir)
    end
    [ny,nx,~] = size(I);
    
    %%%%obtain summed incoherent image
    incoh_img = sum(I(:,:,:,i),3);
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
            contrast = (max_patch-min_patch)/(max_patch+min_patch);
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
            imwrite_float(single(abs(rec_object)),[saveDir 'amplitude.tif']);
            imwrite_float(single(phase),[saveDir 'phase.tif']);
            
        case 'test_ROI'
%             ROIList = [ROI_max];
            ROI(1) = ROI_max(1);
            ROI(2) = ROI_max(3);
            ROI(3) = sx;
            ROI(4) = sx;
            disp(['reconstructing max contrast ROI: [' num2str(ROI(1)) ' ' ...
                 num2str(ROI(2)) ' ' num2str(ROI(3)) ' ' num2str(ROI(4)) ']'])
            [rec_object,phase,rec_pupil,~,~,~,~,~] = AlgorithmManual(...
                I, LEDs, imageColOrder, LEDsUsed, ROI, [], ...
                systemSetup, options, others, tmp0.svdIdx, 2);
            imwrite_float(single(abs(rec_object)),[saveDir 'test_amplitude.tif']);
            imwrite_float(single(phase),[saveDir 'test_phase.tif']);
            
            %I, LEDs, imageColOrder, LEDsUsed, ROIList(90,:), ...
        case 'all_ROIs'
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
    
    
    %%%%do not pass ROI_bg for now, passing 2 to recType means it's not
    %%%%synthetic data
    
end
end

function ROI = centralROI(nx,ny,sx,sy)
    ROI = zeros(4,1);
    ROI(1) = floor(ny/2) - ceil(sy/2);
    ROI(2) = floor(nx/2) - ceil(sx/2);
    ROI(3) = sy;
    ROI(4) = sx;
end
