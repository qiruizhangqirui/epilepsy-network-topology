function dice_martix=image_cross_dice(IG1,IG2,mask,IG2_atlas)
% Function of 'Dice_between_2group' 
% IG1: a cell include all nii/nii.gz path of image group 1
% IG2: a cell include all nii/nii.gz path of image group 2
% mask: nii file of brain mask
% IG2_atlas: 1. Use a network atlas nii/nii.gz path(char) instead of IG2.  Atlas should use 
%               positive integers to identify each brain partition.Background marked with 0.
%            0. Use a group of image for IG2
    info = niftiinfo(mask);
    V = niftiread(info);
    V_2d = reshape(V,[size(V,1)*size(V,2)*size(V,3),1]);
    V_2d_mask = zeros(size(V_2d));
    V_2d_mask(find(V_2d > 0)) = 1;


if IG2_atlas == 1


    info = niftiinfo(IG2);
    V = niftiread(info);
    V_2d = reshape(V,[size(V,1)*size(V,2)*size(V,3),1]);
    atlas_masked_data = V_2d(find(V_2d_mask ==1 ));

    for i = 1:length(IG1)
        
        info2 = niftiinfo(IG1{i});
        V2 = niftiread(info2);
        V2_2d = reshape(V2,[size(V,1)*size(V,2)*size(V,3),1]);

        IG1_data = V2_2d(find(V_2d_mask ==1 ));
        IG1_data(find(IG1_data > 0)) = 1;
        IG1_data(find(IG1_data < 0)) = 0;

        network_label = unique(atlas_masked_data);

        if network_label(1) == 0
            network_label(1) = [];
        end

        for k = 1:length(network_label)

            atlas_data_network = zeros(size(atlas_masked_data));
            atlas_data_network(find(atlas_masked_data == network_label(k))) = 1;
            dice_martix(i,k) = dice(double(IG1_data),double(atlas_data_network));
        end
    end


    else
    for i = 1:length(IG1)
        
        info2 = niftiinfo(IG1{i});
        V2 = niftiread(info2);
        V2_2d = reshape(V2,[size(V,1)*size(V,2)*size(V,3),1]);

        IG1_data = V2_2d(find(V_2d_mask ==1 ));
        IG1_data(find(IG1_data > 0)) = 1;
        IG1_data(find(IG1_data < 0)) = 0;

       
        for k = 1:length(IG2)
        info2 = niftiinfo(IG2{k});
        V2 = niftiread(info2);
        V2_2d = reshape(V2,[size(V,1)*size(V,2)*size(V,3),1]);

        IG2_data = V2_2d(find(V_2d_mask ==1 ));
        IG2_data(find(IG2_data > 0)) = 1;
        IG2_data(find(IG2_data < 0)) = 0;
        dice_martix(i,k) = dice(double(IG1_data),double(IG2_data));
        end
    
    
    
    end
end


