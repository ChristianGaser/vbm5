function p = cg_label_HMRF(p,label,MRFbeta,vx)
% Application of HMRF to segmented data using fast convolution filter
%
% FORMAT p = cg_label_HMRF(p,label,MRFbeta,vx)
% p         - joint probability
% label     - labeled image
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
% @(#)cg_label_HMRF.m	1.06 Christian Gaser 2006/07/14

% default values
if nargin < 3, vx = [1 1 1]; end
if nargin < 2, MRFbeta = 0.3; end
if length(p) ~= 3, error('Number of tissue classes should be 3.'); end

sz = size(p{1});

sum_p = zeros(sz);

for i=1:length(p)
    % compute MRF energy of label
    PUx = get_PUx26(label,i,vx, MRFbeta);

    % calculate joint probabilities P(y,x) = P(y|x)*P(x)
    % P(y|x) are already computed probablities of the tissue classes
    % and P(x) = exp(-U(x)) is the HMRF prior probability of the classes
    p{i} = uint8(round(double(p{i}).*double(PUx)/255));
    sum_p = sum_p + double(p{i});
end

clear PUx;

% Because we don't have the right normalization factor of the Gibbs distribution
% we normalize by forcing that the sum of all joint probabilities is 1
% We have to exclude areas which have low tissue probablity or the overall 
% probability of all tissue classes is too low. 
th = 255*0.05;
mask_label = label > 0;
clear label
for i=1:length(p)
    p{i} = uint8(round(255*double(p{i}).*mask_label.*(p{i}>th)./(sum_p + eps)));
end

return
%=======================================================================

%=======================================================================
function PUx = get_PUx26(label, label_index, vx, MRFbeta)
% Calculation of MRF probability of segmented data
%
% FORMAT PUx = get_PUx26(label, label_index, vx ,MRFbeta)
% label         - labeled segmentation
% label_index   - index of the label to analyze
% vx            - voxel size to correct for anisotropy
% MRFbeta       - bet weighting of MRF energy
%
% Function to calculate prior MRF probability based on MRF energy U(x) according
% to the Ising model by using a convolution filter. The idea consists in using a 
% box of 3x3x3 voxels as filter kernel. Because we are using an image where the 
% voxels of the current label_index are 1 and otherwise 0 this filter calculates 
% the sum of all voxels in the box.
% Distances between the voxel are taken into account by weighting the kernel
% with the invers of the voxel size.
%

% compute an image where voxels with label_index value are 1
sz = size(label);

Ux = repmat(uint8(0),sz);

Ux(find(label==label_index)) = 1;
clear label;

% for anisotropic data we use voxel size to weight kernel according to the distances
% between neighboring voxels
kx = (1/vx(1))*ones(3,1);
ky = (1/vx(2))*ones(3,1);
kz = (1/vx(3))*ones(3,1);

% convolving data with kernel [1 1 1] equals the sum of each label_index
% using 26 neighbours resulting in a maximum value of 9 for all 26 neighbours plus
% the center voxel
spm_conv_vol(Ux,Ux,kx,ky,kz,-[1 1 1]);

% We have to scale the result to make Ux and PUx independent from voxel anisotropy
% For an isotropic cube the maximum is now 9
PUx = double(Ux)*prod(vx)^(2/3);
clear Ux
PUx = exp(MRFbeta*PUx);

% normalize prior probability to 255 (not essential because we
% scale the sum of the joint probability to 1 at the end)
PUx = PUx/max(PUx(:));
PUx = uint8(round(255*PUx));

return
%=======================================================================
