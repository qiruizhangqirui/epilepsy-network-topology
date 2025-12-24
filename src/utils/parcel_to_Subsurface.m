function parcel_to_Subsurface(parcel_value,output_name,color_range,cmap)



f = figure,
    plot_subcortical(parcel_value, 'color_range', color_range, 'cmap', cmap,'ventricles','False')


print(gcf,'-dtiff','-r300',output_name);
close(gcf)

end