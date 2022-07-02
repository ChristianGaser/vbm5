function [g,w,c] = cg_gwc_HMRF(g,w,c,MRFbeta,vx)
% Application of HMRF to segmented data using fast convolution filter
%
% FORMAT [g,w,c] = cg_gwc_HMRF(g,w,c,MRFbeta,vx)
% g         - gray matter image (uint8 format)
% w         - white matter image (uint8 format)
% c         - CSF (uint8 format)
% MRFbeta   - weighting of HMRF
% vx        - voxel size to correct for anisotropy 
%
% Application of a Hidden Markov Random Field (HMRF) model introduces spatial
% constraints based on neighbouring voxels of a 3x3x3 cube. The center voxel has 26
% neighbours and we can calculate MRF energy by counting the number of neighbours.
% The idea is to remove isolated voxels of one tissue class which are unlikely to
% be member of this tissue type. This procedure also closes holes in a cluster
% of connected voxels of one tissue type. In the resulting segmentation the
% noise level will be minimized.
%_______________________________________________________________________
% @(#)cg_gwc_HMRF.m	1.12 Christian Gaser 2006/07/21

fprintf('Applying HMRF with weighting of %3.2f\n', MRFbeta);

sz = size(g);

% use only bounding box where sum of tissue classes is > 0.25 to save memory
th=255*0.25;
mask_gwc = double(g); % separate 'plus' steps to save memory
mask_gwc = mask_gwc + double(w);
mask_gwc = mask_gwc + double(c);
mask_gwc = mask_gwc > th;

[indx, indy, indz] = ind2sub(size(mask_gwc),find(mask_gwc));
indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

% compute MRF energy of label and calculate joint probabilities 
% P(y,x) = P(y|x)*P(x)
% P(y|x) are already computed probablities of the tissue classes
% and P(x) = exp(-U(x)) is the HMRF prior probability of the classes

% GM
label = repmat(uint8(0),sz);
label(find((g >= c) & (g >= w))) = 1;
PUx = get_PUx26(label(indx,indy,indz),vx, MRFbeta);  
g(indx,indy,indz) = uint8(round(double(g(indx,indy,indz)).*double(PUx)/255));
sum_gwc = double(g(indx,indy,indz));

% WM
label = repmat(uint8(0),sz);
label(find((w >= c) & (w >= g))) = 1;
PUx = get_PUx26(label(indx,indy,indz),vx, MRFbeta);
w(indx,indy,indz) = uint8(round(double(w(indx,indy,indz)).*double(PUx)/255));
sum_gwc = sum_gwc + double(w(indx,indy,indz));

% CSF
label = repmat(uint8(0),sz);
label(find((c >= g) & (c >= w))) = 1;
PUx = get_PUx26(label(indx,indy,indz),vx, MRFbeta);
c(indx,indy,indz) = uint8(round(double(c(indx,indy,indz)).*double(PUx)/255));
sum_gwc = sum_gwc + double(c(indx,indy,indz));

clear label PUx;

% Because we don't have the right normalization factor of the Gibbs distribution
% we normalize by forcing that the sum of all joint probabilities is 1
g(indx,indy,indz) = uint8(round(255*double(g(indx,indy,indz)).*mask_gwc(indx,indy,indz)./(sum_gwc + eps)));
w(indx,indy,indz) = uint8(round(255*double(w(indx,indy,indz)).*mask_gwc(indx,indy,indz)./(sum_gwc + eps)));
c(indx,indy,indz) = uint8(round(255*double(c(indx,indy,indz)).*mask_gwc(indx,indy,indz)./(sum_gwc + eps)));

return
%=======================================================================

%=======================================================================
function PUx = get_PUx26(label, vx, MRFbeta)
% Calculation of MRF probability of segmented data
%
% FORMAT PUx = get_PUx26(label, label_index, vx ,MRFbeta)
% label         - labeled segmentation ("1" for current tissue class, "0" for others) 
% vx            - voxel size to correct for anisotropy
% MRFbeta       - beta weighting of MRF energy
%
% Function to calculate prior MRF probability based on MRF energy U(x) according
% to the Ising model by using a convolution filter. The idea consists in using a 
% box of 3x3x3 voxels as filter kernel. Because we are using an image where the 
% voxels of the current tissue class are 1 and otherwise 0 this filter calculates 
% the sum of all voxels in the box.
% Distances between the voxel are taken into account by weighting the kernel
% with the invers of the voxel size.
%

% for anisotropic data we use voxel size to weight kernel according to the distances
% between neighboring voxels
kx = (1/vx(1))*ones(3,1);
ky = (1/vx(2))*ones(3,1);
kz = (1/vx(3))*ones(3,1);

% convolving data with kernel [1 1 1] equals the sum of each label_index
% using 26 neighbours resulting in a maximum value of 9 for all 26 neighbours plus
% the center voxel
spm_conv_vol(label,label,kx,ky,kz,-[1 1 1]);

% We have to scale the result to make Ux (=label) and PUx independent from voxel anisotropy
% For an isotropic cube the maximum is now 9
PUx = double(label)*prod(vx)^(2/3);
clear label
PUx = exp(MRFbeta*PUx);

% normalize prior probability to 255 (not essential because we
% scale the sum of the joint probability to 1 at the end)
PUx = PUx/max(PUx(:));
PUx = uint8(round(255*PUx));

return
%=======================================================================
