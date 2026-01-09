function local_boxplot(groupCat, y, varName, savePath)

    idx = ~isundefined(groupCat) & ~isnan(y);
    figure;
    boxplot(y(idx), groupCat(idx));
    ylabel(varName,'Interpreter','none');
    set(gcf,'Color','w');
    exportgraphics(gcf,savePath,'Resolution',300);

end

