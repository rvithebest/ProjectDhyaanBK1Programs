function set_axis_ticks_fontsize(plotHandles,axis_font_size,tick_font_size,Monkey_num)
    for i=1:size(plotHandles,2)
        subplot(plotHandles(Monkey_num,i))
        ax=gca;
        ax.FontSize=tick_font_size;
        if ~isempty(ax.XLabel.String)
            ax.XLabel.FontSize=axis_font_size;
        end
        if ~isempty(ax.YLabel.String)
            ax.YLabel.FontSize=axis_font_size;
        end
        box(ax, 'off');
    end
end