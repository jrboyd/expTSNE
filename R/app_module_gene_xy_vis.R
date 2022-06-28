server_gene_xy = function(input, output, session, xy_df, color_df, rownames_to_vis, norm_description){
  
  output$plot_tsne_gene = renderPlot({
    req(xy_df())
    req(rownames_to_vis())
    xy = xy_df()
    color_vals = color_df()[rownames_to_vis(),]
    xy$color_val = color_vals[xy$column_id]
    # ggplot(xy, aes(x = tx, y = ty, color = log10(color_val + 1))) +
    ggplot(xy, aes(x = tx, y = ty, color = color_val)) + 
      geom_point() + 
      coord_fixed() +
      labs(x = "", y = "", 
           title = paste(rownames_to_vis(), "expression"), 
           subtitle = norm_description, color = norm_description) +
      scale_color_viridis_c() +
      theme(panel.background = element_blank(),
            panel.grid = element_blank())
    
  })
  
}