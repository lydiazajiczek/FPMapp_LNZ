function AlgorithmManualStitch(ImagesIn, LEDs, imageColOrder, LEDsUsed, ROIList, ROI_bg, saveDir, systemSetup, options, svdIdx,rectype)
% Function that runs FPM algorithms
%   Inputs:
%       ImagesIn - collected images (3d matrix)
%       LEDs - LED matrix used to collect images
%           0 in place where there is no LED
%           1 in place where there is LED
%           2 in place where there is central LED
%           Example:
%               LEDs = [0,1,1,1,0
%                       1,1,2,1,1
%                       0,1,1,1,0];
%       imageColOrder - matrix that shows in which order the images were 
%                       collected
%           1 - first image; 2 - second image; etc
%           Example:
%               imageColOrder = [0,1,2,3,0
%                                4,5,6,7,8
%                               0,9,10,11,0];
%       LEDsUSED - matrix that shows which images you want to reconstruct
%           1 - LED used in reconstruction
%           0 - LED not used in reconstruction
%           Example:
%               LEDsUsed = [0,1,1,1,0
%                           1,1,0,1,1
%                           0,1,1,1,0];
%       ROI - vector of Region Of Interest
%           ROI(i) = [x0,y0,xSize,ySize];
%       systemSetup
%           systemSetup.NA - NA
%           systemSetup.lambda - wavelength (um)
%           systemSetup.magnification magnification
%           systemSetup.LEDspacing - spacing between neighbour LEDs in LED
%                                    matrix (mm)
%           systemSetup.LEDheight - distance between LED matrix and a sample
%                               (mm)
%           systemSetup.camPizSize - camera pixel size (um)
%       options - reconstruction options
%           options.alpha - regularization parameter for object reconstruction
%           options.beta - regularization parameter for pupil reconstruction
%           options.maxIter - maximum number of iterations
%           options.algorithm - select reconstruction algorithm
%               1 - Quasi-Newton algorithm
%               2 - Gerchberg-Saxton
%           options.recorder - reconstruction order
%               1 - from lowest to highest NA
%               2 - from the lightest to the darkest images
%           options.LEDcorrection - LED position correction
%               1 - Angle Self-Calibration
%               2 - simulated annealing
%               3 - genetic algorithm
%           options.initialPupil - initial pupil
%               1 - ones
%               2 - tukey window
%               3 - gauss window
%           options.useGPU
%               1 - use GPU acceleration
%               0 - don't use
%   Outputs:
%       rec_object - reconstructed object (complex double)
%       rec_pupil - reconstructed pupil (complex double)
%       err - RMS error (compared to input data)
%       erro - RMS error (compared to known synthetic object)

%% initialization
tic
% profile on
%f = waitbar(0,'initialization');

%options
initPupil = options.InitPupil;
backgroundROICorr = options.dfBackground;
scaleFactor = options.scaleFactor;
useGPU = options.useGPU;

[ny,nx,nImgs] = size(ImagesIn);
sy = ROIList(1,2)-ROIList(1,1)+1; %size of patch
sx = ROIList(1,4)-ROIList(1,3)+1;

% preparing input images - ROI cropping, background removing, converting to
% GPU array
bck = zeros(nImgs,1);
if ny>100 && nx>100
    for nn = 1:nImgs; bck(nn) = mean2(ImagesIn(1:100,1:100,nn)); end
else
    for nn = 1:nImgs; bck(nn) = mean2(ImagesIn(:,:,nn)); end    
end
thr = (max(bck)+min(bck))/2;

%create stack of background ROIs
if backgroundROICorr
    bgStack = ImagesIn(ROI_bg(2):ROI_bg(2)+ROI_bg(4)-1,ROI_bg(1):ROI_bg(1)+ROI_bg(3)-1,:);
end

%then find initial phase of whole image (to keep phase consistent in stitching)
if options.InitPhase
    disp('generating initial phase estimate (DPC)')
    IDPC = CreateDPCImgs(ImagesIn,LEDsUsed,imageColOrder,false);
    [PhaseIn,~,~] = main_dpc(IDPC,systemSetup,false,options.useGPU);
else
    PhaseIn = zeros(ny,nx);
end
%then crop

%remove background first using ROI of whole image
%if options.dfBackground == 1
%    [ImagesIn, ~] = BackgroundRemovingROI(ImagesIn,ROI_bg,LEDs,imageColOrder,systemSetup);
%else
%    [ImagesIn, ~] = BackgroundRemoving(ImagesIn,thr);
%end

[cledY,cledX] = find(LEDs == 2);    % central LED position
LEDs_old = LEDs; %UGH
LEDs(LEDs>1) = 1;

if LEDsUsed(cledY,cledX) == 0
    cImag = [];
else
    cImag = imageColOrder(cledY,cledX); % central image number
end

%% system setup
lambda = systemSetup.lambda;
NA = systemSetup.NA; 
um_m = NA/lambda;   % maximum spatial frequency set by NA
mag = systemSetup.magnification; 
pixSizeCam = systemSetup.camPizSize; % pixel size on the sensor plane
pixSizeObj = pixSizeCam/mag; % effective image pixel size on the object plane

% Field of view in the object space
FoVx = sx*pixSizeObj;
FoVy = sy*pixSizeObj;

% sampling size in x direction
if mod(sx,2) == 1
    dux = 1/pixSizeObj/(sx-1);
else
    dux = 1/FoVx;
end
% sampling size in y direction
if mod(sy,2) == 1
    duy = 1/pixSizeObj/(sy-1);
else
    duy = 1/FoVy;
end

%% low-pass filter diameter set by the NA = bandwidth of a single measurment
m = 1:sy;
n = 1:sx;
[mm,nn] = meshgrid(m-round((sy+1)/2),n-round((sx+1)/2));
if sy>sx
    nn = nn.*max(max(abs(mm)))./max(max(abs(nn)));
else
    mm = mm.*max(max(abs(nn)))./max(max(abs(mm)));
end
ridx = sqrt(mm.^2+nn.^2);
um_idx = um_m/min(dux,duy);

center = [round(ny/2),round(nx/2)];
[ys,xs] = size(LEDs);
LEDspacing = systemSetup.LEDspacing;    % spacing between adjacent LEDs
LEDheight = systemSetup.LEDheight;  % distance bewteen the LED matrix and the sample
% for LEDheight = 37.5:0.5:42.5
xx = 1:xs; xx = (xx - cledX).*LEDspacing;
yy = 1:ys; yy = (yy - cledY).*LEDspacing;
[LEDsPosX,LEDsPosY] = meshgrid(xx,yy);  % LEDs position in x and y
% LEDsPosX = LEDsPosX - 1;
% distances between LEDs and sample

% % number of brightfield images
% NBF = sum(sum(illuminationNA<NA));
%disp(['synthetic NA is ',num2str(um_p*lambda)]);

%crop in z, not x or y (xy cropping happens inside parfor loop)
%[ImagesIn,imageColOrder] = InputImagesCrop(ImagesIn,imageColOrder,LEDsUsed,[1 1 nx ny]);
%[~, bck] = BackgroundRemoving(ImagesIn,thr); %this is needed for recOrder (intensity)
% reconstruction order
%recOrder = img_order(LEDs, imageColOrder, LEDsUsed, options.recorder, bck);

disp('generating image ROI stacks')
nz = length(ROIList);
I = zeros(sy,sx,nImgs,nz);
PH = zeros(sy,sx,nz);
on_axis = double(ImagesIn(:,:,cImag));
OA = zeros(sy,sx,nz);

%vectorize
for i=1:nz
    coords = ROIList(i,:);
    %colOrder should not change for patches, so overwriting is ok
    [I(:,:,:,i),colOrder] = InputImagesCrop(ImagesIn(coords(1):coords(2),coords(3):coords(4),:),...
                                            imageColOrder,LEDsUsed,[1 1 sx sy]);
    PH(:,:,i) = PhaseIn(coords(1):coords(2),coords(3):coords(4),:);
    OA(:,:,i) = on_axis(coords(1):coords(2),coords(3):coords(4),:);
end

tile_config = zeros(nz,2); %for generating tile config file for IJ stitching

disp('reconstructing patches...')
parfor i=1:nz
%for i=1:2
    %disp('cropping and removing background')
    %[imgs,colOrder] = InputImagesCrop(I(:,:,:,i),imageColOrder,LEDsUsed,[1 1 sx sy]);
    if backgroundROICorr
        [~, bck] = BackgroundRemoving(I(:,:,:,i),thr);
        [I(:,:,:,i), ~] = BackgroundRemovingROI(I(:,:,:,i),bgStack,scaleFactor,LEDs_old,imageColOrder,systemSetup);
    else
        [I(:,:,:,i), bck] = BackgroundRemoving(I(:,:,:,i),thr);
    end
    recOrder = img_order(LEDs, colOrder, LEDsUsed, options.recorder, bck); %#ok<PFBNS>

    coords = ROIList(i,:);
    %I = ImagesIn(coords(1):coords(2),coords(3):coords(4),:); %#ok<PFBNS>
    %PH = PhaseIn(coords(1):coords(2),coords(3):coords(4)); %need to
    %calculate phase from whole image

    %% LEDs position
    img_center_x = (coords(3)-center(1)+sx/2)*pixSizeObj/1000;  %#ok<PFBNS>
    img_center_y = (coords(1)-center(2)+sy/2)*pixSizeObj/1000; %(2)
    
    dist = sqrt((LEDsPosX-img_center_y).^2+(LEDsPosY-img_center_x).^2 ...
    + LEDheight.^2);
    % corresponding angles for each LEDs
    sin_thetaX = (LEDsPosX-img_center_y)./dist;
    sin_thetaY = (LEDsPosY-img_center_x)./dist;

    illuminationNA = sqrt(sin_thetaX.^2+sin_thetaY.^2).*LEDsUsed;
        
    % maxium spatial frequency achievable based on the maximum illumination
    % angle from the LED array and NA of the objective
    um_p = max(max(illuminationNA))/lambda+um_m;
    
    % assume the max spatial freq of the original object
    % um_obj>um_p
    % assume the # of pixels of the original object in x direction
    N_objX2 = round(2*um_p/dux)*2;

    % need to enforce N_obj/Np = integer to ensure no FT artifacts
    N_objX = ceil(N_objX2/sx)*sx;
    if N_objX == sx
        N_objX = N_objX*2;
    end
    N_objY = sy*N_objX/sy;
    tile_config(i,:) = [(coords(3)-1).*uint16(N_objX/sx),(coords(1)-1).*uint16(N_objY/sy)];
    
    % corresponding spatial freq for each LEDs
    xFreq = sin_thetaX/lambda;
    yFreq = sin_thetaY/lambda;
    % spatial freq index for each plane wave relative to the center
    idx_Y = round(yFreq/duy); 
    idx_X = round(xFreq/dux);

    % loading input LEDs positions
    if ~isempty(svdIdx)
        idx_Y = round(svdIdx.idx_Y.*sx./svdIdx.ROI(4));
        idx_X = round(svdIdx.idx_X.*sy./svdIdx.ROI(3));
        idx_Y = idx_Y - idx_Y(cledY,cledX);
        idx_X = idx_X - idx_X(cledY,cledX);
    end
    
    %disp('estimating pupil')
    %then do pupil initialization just of cropped region
    if initPupil == 1
        IDPC = CreateDPCImgs(I(:,:,:,i),LEDsUsed,colOrder,true);
        [~,PupilAmp,PupilPhase] = main_dpc(IDPC,systemSetup,true,options.useGPU);
        phase0 = double(fftshift(PupilPhase));
        pupil0 = double(fftshift(PupilAmp));
        pupil0 = pupil0.*exp(1i.*phase0);
    else
        pupil0 = double(ridx<um_idx);
        pupil0 = initialPupil(pupil0,1);
    end

    if useGPU
        I(:,:,:,i) = gpuArray(I(:,:,:,i));
        %PH(:,:,i) = gpuArray(PH(:,:,i)); %is this needed? prob not
    end 

    % figure; imagesc(pupil0); title ones
    % figure; imagesc(pupil0);  title tukey
    
    %disp('reconstructing')
    %% reconstruction algorithm
    [rec_object, phase, pupil, ~, ~, ~, ~] = ... 
        AlgorithmStitch(I(:,:,:,i), PH(:,:,i), [N_objY,N_objX], idx_X, idx_Y, colOrder, ...
        recOrder, pupil0, options, cImag, rectype); %showIm=false
    %AlgorithmStitch(I, PH, N_obj, idx_X, idx_Y, imageColOrder, recOrder, pupil0, options, cImag, rectype)
    
    %disp('writing to disk')
    %% convert rec_object to uint16 and write to disk
    %% normalize brightness using histogram equalization of on-axis image
    amplitude = abs(rec_object); %convert to real
        
    %normalize to incoherent image (summed and normalized image stack)
    incoh_img = sum(I(:,:,:,i),3); 
    incoh_img = incoh_img.*(max(max(OA(:,:,i)))/max(max(incoh_img)));
    amplitude = uint16((min(min(incoh_img))/(min(min(amplitude))+1)).*amplitude);
    
    %histogram equalization
    [n,~] = histcounts(uint16(incoh_img));
    amplitude = uint16(histeq(amplitude,n));
    imwrite(amplitude,[saveDir '\amplitude' num2str(i,'%03d') '.tif'])
    
    %% normalize phase and write to disk
    phase_in = (PH(:,:,i)+pi)./(2*pi); %convert to [0 1] for histogram normalization
    phase_norm = (phase+pi)./(2*pi);
    [n,~] = histcounts(phase_in);
    phase_norm = (2*pi).*histeq(phase_norm,n); %convert back to [0 2pi]
    
    %my ropey attempt to normalize phase
    %phase_mod = phase.*((max(max(phase_in)) - min(min(phase_in)))/(abs(max(max(phase)) - min(min(phase))))) ...
    %                   - (max(max(phase_in)) - min(min(phase_in)))/2;
    imwrite_float(single(phase),[saveDir '\phase' num2str(i,'%03d') '.tif'])
    
    %save pupil just for fun
    imwrite_float(single(angle(pupil)),[saveDir '\pupil' num2str(i,'%03d') '.tif'])
    
end

%write tile configuration files for ImageJ stitching here
f_amp = fopen([saveDir '\AmplitudeTileConfiguration.txt'],'w+');
fprintf(f_amp,'dim = 2\n');
f_ph = fopen([saveDir '\PhaseTileConfiguration.txt'],'w+');
fprintf(f_ph,'dim = 2\n');
for i=1:nz
   str = ['amplitude' num2str(i,'%03d') '.tif; ; (' num2str(tile_config(:,1)) ',' num2str(tile_config(:,2)) ')\n'];
   fprintf(f_amp,str);
   str = ['phase' num2str(i,'%03d') '.tif; ; (' num2str(tile_config(:,1)) ',' num2str(tile_config(:,2)) ')\n'];
   fprintf(f_ph,str);
end
fclose(f_amp);
fclose(f_ph);

toc
end

