using CairoMakie

my_thesis_theme = Theme(fontsize=12,
                        font="TeX Gyre Heros",
                        Lines=(linewidth=2,),
                        palette=(color=to_colormap(:seaborn_bright6),),
                        Axis=(xlabelsize=14,
                              ylabelsize=14,
                              titlesize=16,
                              xticklabelsize=12,
                              yticklabelsize=12),
                        Legend=(fontsize=12,))

cm_to_px(cm) = cm / 2.54 * 72
TEXT_WIDTH_CM = 16.5
TEXT_WIDTH_PX = cm_to_px(TEXT_WIDTH_CM)

function my_figure_size(width_px; aspect_ratio=1.618)
    return (width_px, width_px / aspect_ratio)
end

function export_fig(file_name, fig; png_dpi=300, save_pdf=false, save_png=false)
    base_dpi = 72
    dir = dirname(file_name)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end
    if save_pdf
        save("$(file_name).pdf", fig)  # Vektor – skaliert verlustfrei
    end
    if save_png
        scale = png_dpi / base_dpi     # z.B. 300/72 ≈ 4.17
        save("$(file_name).png", fig; px_per_unit=scale)
    end
end
