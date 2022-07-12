

ui_point_selection = function(id = "default_point_selection", ps_size = "400px", n_items = 10){
  ns = NS(id)
  tabsetPanel(
    id = ns("tabs_cluster_or_selection"),
    tabPanel(
      "Clustering",
      sidebarLayout(
        sidebarPanel(
          tags$h5("TODO: clean up some weird A/B interaction, B trumps A currently."),
          tags$label("Patient Clustering Methods"),
          tabsetPanel(
            id = ns("tabs_cluster_method"),
            tabPanel(
              "knn",
              numericInput(ns("num_nn"), label = "Nearest neighbors", value = ceiling(n_items/5), min = 2, max = Inf)
            ), 
            tabPanel(
              "kmeans",
              numericInput(ns("num_kmeans"), label = "k", value = min(5, n_items-1), min = 2, max = Inf)
            ), 
            tabPanel(
              "hierarchical",
              numericInput(ns("num_clust"), label = "n_clust", value = min(5, n_items-1), min = 2, max = Inf)
            )
          )), 
        mainPanel(
          shinycssloaders::withSpinner(plotOutput(ns("plot_tsne_clusters"), width = ps_size, height = ps_size))  
        )
      )),
    tabPanel(
      "Selection",
      sidebarLayout(
        sidebarPanel(
          shinyjs::disabled(selectInput(ns("sel_A_clust"), label = "A group clusters", choices = "", multiple = TRUE)),
          shinyjs::disabled(selectInput(ns("sel_B_clust"), label = "B group clusters", choices = "", multiple = TRUE)),    
        ),
        mainPanel(
          shinycssloaders::withSpinner(plotOutput(ns("plot_A_group"), width = ps_size, height = ps_size,
                                                  brush = brushOpts(id = ns("plot_A_brush")))),
          tags$h3("Refine selection"),
          fluidRow(
            column(width = 4,
                   actionButton(ns("btn_add_A"), "Add A"),
                   actionButton(ns("btn_rm_A"), "Remove A"),
                   actionButton(ns("btn_limit_A"), "Limit A"),
            ), 
            column(width = 4,
                   actionButton(ns("btn_add_B"), "Add B"),
                   actionButton(ns("btn_rm_B"), "Remove B"),
                   actionButton(ns("btn_limit_B"), "Limit B")
            )
            
          ),
        )
      )
    ),
  )
}

server_point_selection = function(id = "default_point_selection", 
                                  tsne_clust, 
                                  meta_data, 
                                  tsne_res, 
                                  column_id_var = "column_id",
                                  x_var = "tx",
                                  y_var = "ty",
                                  cluster_var = "tsne_cluster_id"){
  moduleServer(
    id,
    function(input, output, session){
      output$plot_tsne_clusters = renderPlot({
        req(tsne_clust())
        tsne_dt = tsne_clust()
        tsne_dt = merge(tsne_dt, meta_data(), by = column_id_var)
        
        p = ggplot(tsne_dt, aes_string(x = x_var, y = y_var, color = cluster_var)) + 
          geom_point() + 
          coord_fixed() +
          labs(x = "", y = "", title = "patient clustering", subtitle = paste(input$sel_gene_list, "gene list")) +
          theme(panel.background = element_blank(),
                panel.grid = element_blank())
        p
      })
      
      observe({
        if(is.null(tsne_res())){
          tsne_clust(NULL)
        }
        req(tsne_res())
        message("run clustering")
        tsne_dt = as.data.table(tsne_res())
        clust_method = input$tabs_cluster_method
        showNotification(paste("clustering method is", clust_method))
        if(clust_method == "knn"){
          tsne_dt.clust = nn_clust(tsne_dt, nn = input$num_nn, cluster_var = cluster_var)
        }else if(clust_method == "kmeans"){
          tsne_dt.clust = km_clust(tsne_dt, k = input$num_kmeans, cluster_var = cluster_var)
        }else if(clust_method == "hierarchical"){
          tsne_dt.clust = h_clust(tsne_dt, n_clust = input$num_clust, cluster_var = cluster_var)
        }else{
          stop("Unrecognized clustering method, ", clust_method)
        }
        tsne_dt.clust$group = "bg"
        tsne_clust(tsne_dt.clust)
        # updateSelectInput(session, "sel_A_clust", choices = cl, selected = cl[1])
        # updateSelectInput(session, "sel_B_clust", choices = cl, selected = cl[2])
      })
      
      observeEvent({
        tsne_clust()
      },
      {
        if(is.null(tsne_clust())){
          disable("sel_A_clust")
          disable("sel_B_clust")
        }
        cl = sort(unique(tsne_clust()[[cluster_var]]))
        cl = cl[order(as.numeric(sub("[a-zA-Z ]+", "", cl)))]
        shinyjs::enable("sel_A_clust")
        shinyjs::enable("sel_B_clust")
        updateSelectInput(session, "sel_A_clust", choices = cl, selected = cl[1])
        updateSelectInput(session, "sel_B_clust", choices = cl, selected = cl[2])
        prev_A_clust = cl[1]
        prev_B_clust = cl[2]
        sel_A_ids(tsne_clust()[get(cluster_var) %in% cl[1]][[column_id_var]])
        sel_B_ids(tsne_clust()[get(cluster_var) %in% cl[2]][[column_id_var]])
        
      })
      
      
      sel_A_ids = reactiveVal(character())
      sel_B_ids = reactiveVal(character())
      prev_A_clust = character()
      prev_B_clust = character()
      
      update_ids = function(current_clust, previous_clust, current_ids){
        added_clust = setdiff(current_clust, previous_clust)
        removed_clust = setdiff(previous_clust, current_clust)
        if(length(added_clust) > 0 & length(removed_clust) > 0){
          stop("inconceivable! new and missing longer than 0.")
        }
        tsne_dt = tsne_clust()
        new_ids = current_ids
        if(length(added_clust) > 0){
          new_ids = union(
            current_ids,
            tsne_dt[get(cluster_var) %in% added_clust][[column_id_var]]
          )
        }
        if(length(removed_clust) > 0){
          new_ids = setdiff(
            current_ids,
            tsne_dt[get(cluster_var) %in% removed_clust][[column_id_var]]          
          )
        }
        new_ids
      }
      
      observeEvent({
        input$sel_A_clust
      }, {
        new_ids = update_ids(input$sel_A_clust, 
                             prev_A_clust,
                             isolate(sel_A_ids()))
        prev_A_clust = input$sel_A_clust
        sel_A_ids(new_ids)
      })
      
      observeEvent({
        input$sel_B_clust
      }, {
        new_ids = update_ids(input$sel_B_clust, 
                             prev_B_clust,
                             isolate(sel_B_ids()))
        prev_B_clust = input$sel_B_clust
        sel_B_ids(new_ids)
      })
      
      output$plot_A_group = renderPlot({
        tsne_dt = tsne_clust()
        tsne_dt$group = "."
        tsne_dt[get(column_id_var) %in% sel_A_ids(), group := "A"]
        tsne_dt[get(column_id_var) %in% sel_B_ids(), group := "B"]
        
        tsne_dt$group = factor(tsne_dt$group, )
        
        p = ggplot(tsne_dt, aes_string(x = x_var, y = y_var, color = "group")) +
          geom_point(data= tsne_dt[group == "."], size = .3) +
          geom_point(data= tsne_dt[group != "."], size = .8) +
          scale_color_manual(values = c("A" = "red", "B" = "blue", "." = "gray"), drop = FALSE) +
          coord_fixed() +
          theme(panel.background = element_blank(),
                panel.grid = element_blank()) +
          labs(x = "", y = "", title = "A and B selection")
        p
      })
      
      filter_by_brush = function(tsne_dt, brsh, column_id_var){
        tsne_dt[get(x_var) >= brsh$xmin & 
                  get(x_var) <= brsh$xmax & 
                  get(y_var) >= brsh$ymin & 
                  get(y_var) <= brsh$ymax][[column_id_var]]
      }
      
      #plot A buttons
      observeEvent({
        input$btn_add_A
      }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = filter_by_brush(tsne_dt, brsh, column_id_var)
        
        sel_A_ids(
          union(isolate(sel_A_ids()), 
                ids_in_rng)    
        )
      })
      
      
      
      observeEvent({
        input$btn_rm_A
      }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = filter_by_brush(tsne_dt, brsh, column_id_var)
        
        sel_A_ids(
          setdiff(isolate(sel_A_ids()), 
                  ids_in_rng)    
        )
      })
      
      observeEvent({
        input$btn_limit_A
      }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = filter_by_brush(tsne_dt, brsh, column_id_var)
        
        sel_A_ids(
          intersect(isolate(sel_A_ids()), 
                    ids_in_rng)    
        )
      })
      
      observeEvent({
        input$btn_add_B
      }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = filter_by_brush(tsne_dt, brsh, column_id_var)
        
        sel_B_ids(
          union(isolate(sel_B_ids()), 
                ids_in_rng)    
        )
      })
      
      observeEvent({
        input$btn_rm_B
      }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = filter_by_brush(tsne_dt, brsh, column_id_var)
        
        sel_B_ids(
          setdiff(isolate(sel_B_ids()), 
                  ids_in_rng)    
        )
      })
      
      observeEvent({
        input$btn_limit_B
      }, {
        brsh = input$plot_A_brush
        tsne_dt = tsne_clust()
        ids_in_rng = filter_by_brush(tsne_dt, brsh, column_id_var)
        sel_B_ids(
          intersect(isolate(sel_B_ids()), 
                    ids_in_rng)    
        )
      })
      
      # id_groups = reactiveValues(A = list(), B = list())
      
      A = reactiveVal()
      B = reactiveVal()
      
      observeEvent({
        sel_A_ids()
        sel_B_ids()
      }, {
        A(sel_A_ids)
        B(sel_B_ids)
      })
      
      list(A = A, B = B)
    })
}