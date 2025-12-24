function parcel_to_volume(parcel_value,atlas,output_name,atlas_order)


    info = niftiinfo(atlas);
    V = niftiread(info);
    V_2d = reshape(V,[size(V,1)*size(V,2)*size(V,3),1]);
    V1_2d = zeros(size(V_2d));
    parcel_index = unique(V(:));
    for i = 1:length(parcel_index)-1
     ind = find(V_2d == atlas_order(i));
     V1_2d(ind) = parcel_value(i);
    end
    V1_3d = reshape(V1_2d,[size(V,1),size(V,2),size(V,3)]);
    if contains(info.Datatype,'int16' )    
        V1_3d = int16(V1_3d);
    elseif contains(info.Datatype,'single' )    
        V1_3d = single(V1_3d);
           elseif contains(info.Datatype,'uint8' )    
        info.Datatype = 'double'
    end
%     info= rmfield(info, 'Intent')
% info= rmfield(info, 'IntentDescription')
% info= rmfield(info, 'IntentParams')
    niftiwrite(V1_3d, output_name, info);
    end