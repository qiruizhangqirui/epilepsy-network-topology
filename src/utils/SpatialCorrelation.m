function [r,p,r_dist,p_spin,Image1,Image2]=SpaitalCorrelation(IM1,IM2,BrainMask,nPer)
% BrainMask = mask;
% [r,p,r_dist,p_spin]=SpaitalCorrelation(IM1,IM2,BrainMask,nPer)
%   
% Inputs£º
%   IM1,IM2        the directory of Nifti images or N * 1 matrix
%   BrainMask     the directory of BrainMask or  N * 1 matrix
%   nPer              times of spin permutation test,0 indicate without spin permutation test
%                           
% Outputs:
%   r                      r value of pearson corrrealtion analysis
%   p                     p value of pearson corrrealtion analysis
%   r_dist              distribution of r values for spin permutation test
%   p_spin             p value of spin permutation test
%
%  References:
%   Alexander-Bloch A, Shou H, Liu S, Satterthwaite TD, Glahn DC, Shinohara RT,
%       Vandekar SN and Raznahan A (2018). On testing for spatial correspondence
%       between maps of human brain structure and function. NeuroImage, 178:540-51.
%   Larivi¨¨re, S., Paquola, C., Park, By. Royer, J., Wang, Y., Benkarim, O., Vos de Wael, R., 
%       Valk, S., Thomopoulos, S.I., Kirschner, M., Lewis, L.B., Evans, A.C., Sisodiya, S.M., 
%       McDonald, C.R., Thompson, P.T, Bernhardt, B.C.. The ENIGMA Toolbox: multiscale
%       neural contextualization of multisite neuroimaging datasets. Nat Methods 18, 698¨C700 

if ischar(IM1) || ischar(IM2)
    I1 = y_Read(IM1);Image1 = I1(:);
    I2 = y_Read(IM2);Image2 = I2(:);
elseif size(IM1, 1) == 1
    Image1 = IM1.';
    Image2 = IM2.';
elseif nnz(size(IM1)) >2 || nnz(size(IM2)) >2 
    disp('The matrix size  is not N*1, please check it !!')
else
    Image1 = IM1;
    Image2 = IM2;
end

if isempty(BrainMask)
    disp('Do not apply brainmask!')
elseif isnumeric(BrainMask)
    Image1 = Image1(BrainMask~=0);Image2 = Image2(BrainMask~=0);
elseif ischar(BrainMask)
    BrainMask = y_Read(BrainMask);
    BrainMask = BrainMask(:);
    BrainMask(Image1 == 0) =0;
    BrainMask(Image2 == 0) =0;
    Image1 = Image1(BrainMask ~=0);Image2 = Image2(BrainMask ~=0);
end

[r,p] = corr(Image1,Image2,'type','Pearson','rows','pairwise');

if nPer >0
    [r_dist,p_spin] = Spin_test(Image1,Image2,nPer);
end

return

function [r_dist,p_spin] = Spin_test(Image1,Image2,nPer)
% build up the matrix of permutation index
perm_id = zeros(length(Image1),nPer);

Times=0;Overlay=0;
while ( Times < nPer)
    
    Temp1 = randperm(length(Image1))';
    
    % get the permutation index
    if ~all(Temp1'==( 1:length(Image1) ))
        Times = Times+1;
        perm_id(:,Times) = Temp1;
    else % verify that permutation does not map to itself
        Overlay = Overlay+1;
        disp(['map to itself n.' num2str(Overlay)])
    end
    
    % Track progress
    if mod(Times,100)==0
        disp(['permutation ' num2str(Times) ' of ' num2str(nPer)]);
    end
    
end

% empirical correlation
rho_emp = corr(Image1, Image2, 'type','Pearson', 'rows', 'pairwise');

% permutation of images
for i =1:nPer
    x_perm(:,i) = Image1(perm_id(:,i));
    y_perm(:,i) = Image2(perm_id(:,i));
end

% corrrelation to unpermuted measures
rho_null_xy = zeros(nPer, 1);
rho_null_yx = zeros(nPer, 1);
for r = 1:nPer
    rho_null_xy(r,1) = corr(x_perm(:, r), Image2, 'type', 'Pearson', 'rows', 'pairwise'); % correlate permuted x to unpermuted y
    rho_null_yx(r,1) = corr(y_perm(:, r), Image1, 'type', 'Pearson', 'rows', 'pairwise'); % correlate permuted y to unpermuted x
end

% p-value definition depends on the sign of the empirical correlation
if rho_emp > 0
    p_perm_xy = sum(rho_null_xy > rho_emp)/nPer;
    p_perm_yx = sum(rho_null_yx > rho_emp)/nPer;
else
    p_perm_xy = sum(rho_null_xy < rho_emp)/nPer;
    p_perm_yx = sum(rho_null_yx < rho_emp)/nPer;
end
% null distribution
r_dist = [rho_null_xy; rho_null_yx];

% average p-values
p_spin = (p_perm_xy + p_perm_yx)/2;
return