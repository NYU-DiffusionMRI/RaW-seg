% Random Walker segmentation (RaW)
% author: Hong-Hsi Lee, 2018

clear all
root = '/data1/Hamster/Honghsi/Projects/segmentation3D/share_code';
rootdata = fullfile(root,'data');
addpath(genpath(fullfile(root,'tools')));
addpath(genpath(fullfile(root,'lib')));
target = fullfile(root,'result'); mkdir(target);

%% Calculate foreground mask and correct mild distortion
% You need to download our EM images to run this section.

target_proc = fullfile(root,'processed'); mkdir(target_proc)
data = load('/data1/Hamster/Honghsi/Projects/segmentation3D/allData.mat');
data = data.data(:,:,1:200);

% Step 1. Use magicwand to create foreground mask
% Input:
%   data: EM data
%   seed: seeding positions. It is in the same size of the data. You need
%     to apply a few seeds in the background for each slice.
rs = rawseg();
seed = false(size(data));
seed(:,1,:) = true;  % apply a few seeds in the background
mask = rs.magicwand(data,seed);
save(fullfile(target_proc,'foregroundmask.mat'),'mask')

% Step 2. Use optical flow to correct mild distortion
% Input: EM data and foreground mask
% Output:
%   datac: corrected data
%   I: the slice with distortion
%   ccf: correlation coefficient btw distorted slice and its interpolation
%   maskc: corrected mask
[datac,Ic,ccf,maskc] = rs.distortioncorrect(data,mask);
save(fullfile(target_proc,'dataDistortionCorrect.mat'),'datac','Ic','ccf','maskc')
save(fullfile(rootdata,'foregroundmaskDistortionCorrect.mat'),'maskc');

% Step 3. Use pixel-wise classifier (e.g., ilastik) to create the myelin
% mask. You may need to save datac into hdf5 format before using the
% classifier.

%% RaW Segmentation

% Load myelin mask
maskmy = load(fullfile(rootdata,'myelinmask.mat'));
maskmy = maskmy.mask;

% Load background mask
maskfg = load(fullfile(rootdata,'foregroundmaskDistortionCorrect.mat'));
maskbg = ~maskfg.maskc;
% Dilate backgound mask to avoid edge effect
maskbg = imdilate(maskbg,strel('cuboid',[15,15,1]));

% Combine two masks to create medium for RaW segmentation
medium = logical(maskmy + maskbg);
save(fullfile(target,'medium.mat'),'medium')

% Load seeding positions. You have to do it manually.
seed = load(fullfile(rootdata,'seed.mat'));
seed = seed.seed;

% Random hopping from seeding positions, one seed per fiber
rs = rawseg();
Npar = 4e3;  % particle number
Nstep = 16*size(medium,3).^2;  % step number
mkdir(fullfile(target,'step1_randomhopping'));
tic;
parfor i = 1:size(seed,1)
    fiber = rs.randomhopping(medium,seed(i,:),Npar,Nstep);
    fiber = logical(fiber);
    rs.savefiber(fullfile(target,'step1_randomhopping',sprintf('fiber%u.mat',i)),fiber);
end
toc;

%% The first proofreading
% Check the file size: The fiber with leaky myelin mask has larger volume,
% which is proportional to the file size.
files = dir(fullfile(target,'step1_randomhopping','fiber*.mat'));
C = struct2cell(files);
filesize = cell2mat(C(4,:));
filesize = filesize(:)/2^20; % file size, MB
figure; hist(filesize,100);  % histogram shows that files < 3 MB are good
xlabel('file size (MB)'); ylabel('frequency')
I = [];
for i = 1:numel(files)
    if filesize(i) < 3
        filename = files(i).name;
        Ii = str2double(filename(6:end-4));
        I = cat(1,I,Ii);
    end
end
I = sort(I);
save(fullfile(target,'proofread1_fiberlabel.mat'),'I')

%% random hopping from seeds chosen from the previous segmentation
load(fullfile(target,'proofread1_fiberlabel.mat'))
load(fullfile(target,'medium.mat'))
rs = rawseg(); Npar = 1e4; Nstep = 4e4;
mkdir(fullfile(target,'step2_randomhopping'));
tic;
parfor i = 1:numel(I)
    filename = sprintf('fiber%u.mat',I(i));
    fiber = load(fullfile(target,'step1_randomhopping',filename));
    fiber = fiber.fiber;
    [nx,ny,~] = size(fiber);
    zlist = find(sum(sum(fiber,1),2));
    % choose seeds
    nseed = numel([zlist(1:10:end).',zlist(end)]);
    seed = zeros(nseed,3);
    j = 0;
    for k = [zlist(1:10:end).',zlist(end)]
        j = j+1;
        [x,y] = ind2sub([nx,ny],find(squeeze(fiber(:,:,k))));
        [~,Ik] = datasample(x,1);
        seed(j,:) = [x(Ik),y(Ik),k];
    end
    % random hopping
    fiber = rs.randomhopping(medium,seed,Npar,Nstep);
    fiber = logical(fiber);
    rs.savefiber(fullfile(target,'step2_randomhopping',filename),fiber);
end
toc;

%% The second proofreading
mkdir(fullfile(target,'proofread2'))
load(fullfile(target,'proofread1_fiberlabel.mat'))
nfig = ceil(numel(I)/20);
rs = rawseg();
for i = 1:nfig
    close all
    h = figure;
    j = 0;
    for k = (i-1)*20+1:min(i*20,numel(I))
        j = j+1;
        filename = sprintf('fiber%u.mat',I(k));
        fiber = load(fullfile(target,'step2_randomhopping',filename));
        subplot(4,5,j)
        rs.visualizefiber(fiber.fiber);
        title(sprintf('fiber%u',I(k)));
    end
    savefig(h,fullfile(target,'proofread2',sprintf('group%u.fig',i)))
end

% After proofreading the fiber shape, we noticed that fiber 1 = fiber 7,
% fiber 2 = fiber 9, and fiber 292 and 322 are bizarre.
I = setdiff(I,[7,9,292,322]);
save(fullfile(target,'proofread2_fiberlabel.mat'),'I')

% Save all fiber segmentations into one file
load(fullfile(target,'proofread2_fiberlabel.mat'))
fiber = load(fullfile(target,'step2_randomhopping',sprintf('fiber%u.mat',I(1))));
fibers = zeros(size(fiber.fiber),'uint16');
masks = zeros(size(fiber.fiber),'uint8');
parfor i = 1:numel(I)
    filename = sprintf('fiber%u.mat',I(i));
    fiber = load(fullfile(target,'step2_randomhopping',filename));
    fibers = fibers + uint16(I(i)*fiber.fiber);
    masks = masks + uint8(fiber.fiber>0);
end
save(fullfile(target,'fibers.mat'),'fibers')
save(fullfile(target,'masks.mat'),'masks')

%% Check the final segmentation by combining it with the EM image.
% You may need to save fibers into hdf5 format before using other software 
% (e.g., ilastik).
% Here, we create a video of the segmentation overlaid on the EM image.
load(fullfile(root,'processed','dataDistortionCorrect.mat'));
load(fullfile(target,'fibers.mat'))
rs = rawseg();
rs.animatefiber(fullfile(target,'fibers.avi'),datac,fibers);

% We also confirmed that no segmentations overlap with each other.
load(fullfile(target,'masks.mat'))
if max(masks(:)) > 1
    fprintf('Segmentations overlap with each other.\n')
else
    fprintf('No segmentations overlap with each otehr.\n')
end

%% Resize matrix to isotropic resolution and fill the holes in each slice
load(fullfile(target,'fibers.mat'))
load(fullfile(target,'proofread2_fiberlabel.mat'))

vox = [24,24,100]*1e-3;  % voxel size, µm
rs = rawseg();
fiberiso = rs.resize(fibers,vox);

fiberfill = zeros(size(fiberiso),'uint16');
parfor i = 1:numel(I)
    fiberi = fiberiso == I(i);
    fiberi = rs.fillhole(fiberi);
    fiberfill = fiberfill + uint16(fiberi*I(i));
end
save(fullfile(target,'fiberfill.mat'),'fiberfill')

%% Segment myelin sheath for each fiber
load(fullfile(target,'fiberfill.mat'))

% Watershed
rs = rawseg();
L = rs.watershed(fiberfill);
save(fullfile(target,'watershed.mat'),'L')
%%
% Dilate IAS and segment individual myelin sheath
vox = [24,24,100]*1e-3;  % voxel size, µm
rs = rawseg();
maskmy = load(fullfile(rootdata,'myelinmask.mat')); 
maskmy = rs.resize(single(maskmy.mask),vox) > 0.5;
maskfg = load(fullfile(rootdata,'foregroundmaskDistortionCorrect.mat'));
maskfg = rs.resize(single(maskfg.maskc),vox) > 0.5;

load(fullfile(target,'fiberfill.mat'))
load(fullfile(target,'proofread2_fiberlabel.mat'))
load(fullfile(target,'watershed.mat'))

vox = [100,100,100]*1e-3;  % voxel size, µm
myelinmax = 0.4;           % maximal myelin thickness, µm

rs = rawseg();
myelins = rs.myelinsheath(fiberfill,I,L,maskmy,maskfg,myelinmax,vox);
save(fullfile(target,'myelinsheath.mat'),'myelins');

%%
addpath(genpath('/data1/Hamster/Honghsi/Projects/segmentation3D/tools/NIfTI_tool'))
nii = make_nii(uint16(fibers),[24,24,100]*1e-3,[],4);
save_nii(nii,fullfile(root,'processed','fibers.nii'));


