ui_goi = function(id = "default_goi"){
  ns = NS(id)
  sidebarLayout(
    sidebarPanel(
      uiOutput(ns("ui_sel_color_by")),
      selectInput(ns("txtGene"), label = "Select Gene To Plot", choices = NULL)
    ),
    mainPanel(
      withSpinner(plotOutput(ns("plot_tsne_gene"), width = "600px", height = "600px"))
    )
  )
}

server_goi = function(id = "default_goi", expression_data, xy_df){
  moduleServer(
    id, 
    function(input, output, session){
      vis_gene = reactiveVal()
      observeEvent({
        expression_data()
      }, {
        req(expression_data())
        def = vis_gene()
        if(is.null(def)) def = ""
        all_genes = rownames(expression_data())
        if(def %in% all_genes){
          updateSelectizeInput(session, 'txtGene', choices = all_genes, selected = def, server = TRUE)
        }else{
          updateSelectizeInput(session, 'txtGene', choices = all_genes, server = TRUE)
        }
        
      })
      observeEvent({
        input$txtGene   
        expression_data()
      }, {
        req(expression_data())
        if(input$txtGene %in% rownames(expression_data())){
          vis_gene(input$txtGene)
        }
      })
      
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
      
    })
}