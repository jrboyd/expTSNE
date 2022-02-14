server_gene_xy = function(input, output, session, xy_df, meta_data, rownames_to_vis){
    
    output$plot_tsne_gene = renderPlot({
        req(xy_df())
        xy = xy_df()
        req(rownames_to_vis())
        # browser()
        xy = merge(xy, meta_data(), by = "column_id")
        color_val = "group"
        # color_vals = meta_data()[rownames_to_vis(),]
        # xy$color_val = color_vals[xy$patient_id]
        # ggplot(xy, aes(x = x, y = y, color = log10(color_val + 1))) + 
        #     geom_point() + 
        #     coord_fixed() +
        #     labs(x = "", y = "", title = paste(rownames_to_vis(), "expression"), subtitle = "log10 scale") +
        #     scale_color_viridis_c() +
        #     theme(panel.background = element_blank(),
        #           panel.grid = element_blank())
        ggplot(xy, aes_string(x = "tx", y = "ty", color = color_val)) + 
            geom_point() + 
            coord_fixed() +
            labs(x = "", y = "", subtitle = "log10 scale") +
            # scale_color_viridis_c() +
            theme(panel.background = element_blank(),
                  panel.grid = element_blank())
        
    })
    
}