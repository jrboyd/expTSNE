
expTSNE.ui_module = function(id){
  ns = NS(id)
  #ui wrap in ns()
  tags$span(
    shinyjs::useShinyjs(),
    tabsetPanel(
      tabPanel(
        "Main",
        # Sidebar with a slider input for number of bins 
        sidebarLayout(
          sidebarPanel(
            uiOutput(ns("ui_sel_sample_type_filter")),
            uiOutput(ns("ui_sel_facet_var")),
            uiOutput(ns("ui_sel_gene_list")),
            disabled((selectInput(ns("sel_custom_gene_set"), label = "Custom gene lists", choices = ""))),
            uiOutput(ns("ui_sel_color_by")),
            selectInput(ns("txtGene"), label = "Select Gene To Plot", choices = NULL),
          ),
          mainPanel(
            withSpinner(plotOutput(ns("plot_tsne"), width = "600px", height = "600px")),
            withSpinner(plotOutput(ns("plot_tsne_gene"), width = "600px", height = "600px"))
          )
        )
      ),
      tabPanel(
        "Add Gene Set",
        ui_upload(id = ns("default_upload"))
      ),
      tabPanel(
        "DE",
        ui_point_selection(id = ns("default_point_selection")),
        tabsetPanel(
          tabPanel("simple",
                   actionButton(ns("btn_runDEfast"), "Run DE"),
                   withSpinner(plotOutput(ns("plot_volcano"), width = "400px", height = "400px")),
                   withSpinner(plotOutput(ns("plot_boxes"), width = "400px", height = "400px"))
                   
          ),
          tabPanel("DESeq2",
                   actionButton(ns("btn_runDE"), "Run DE"),
                   withSpinner(DT::dataTableOutput(ns("dt_DE_res")))
          )
        )
      )
    )
  )
}

expTSNE.server_module = function(id, et){
  moduleServer(
    id,
    function(input, output, session) {
      #server logic
      FACET_VAR = list(NONE = "none", SAMPLE_TYPE = "sample type", PAM50 = "PAM50")
      gene_lists = list(
        All = rownames(et$norm_counts)
      )
      ### TCGA expression data
      expression_data = reactiveVal(et$norm_counts)
      tsne_input = reactiveVal(et$norm_counts)
      ### Genes to use in t-sne
      input_genes = reactiveVal(rownames(et$norm_counts))
      #subset of parsed genes in tcga expression data
      valid_genes = reactiveVal()
      #table shown in Add Gene Set main panel, includes parsed genes and other data used for selection
      gene_table = reactiveVal(data.frame(gene_name = rownames(et$norm_counts)))
      #all added custom gene sets
      custom_gene_sets = reactiveVal(list())
      
      #metadata for patient entries
      meta_data = reactiveVal(et$meta_data)
      #metadata for sample entries
      active_dataset = reactiveVal(et)
      
      #result of tsne
      tsne_res = reactiveVal(as.data.table(et$tsne_result))
      #tsne_res with clustering applied
      tsne_clust = reactiveVal()
      
      dataset_downstream = list(
        # sample_groups_A, 
        # sample_groups_B, 
        # DE_res, 
        # DE_fast_raw, 
        # DE_fast_res,
        # DE_res
      )
      
      output$ui_sel_facet_var = renderUI({
        ns <- session$ns
        dt = meta_data()
        req(dt)
        n_unique = sapply(colnames(dt), function(nam){
          length(unique(dt[[nam]]))
        })
        facet_vars = names(n_unique)[n_unique <= 16]
        radioButtons(
          inputId = ns("sel_facet_var"), 
          label = "Facet By", 
          choices = facet_vars)
      })
      
      output$ui_sel_gene_list = renderUI({
        ns <- session$ns
        tagList(
          radioButtons(
            inputId = ns("sel_gene_list"), 
            label = "Gene List", 
            choices = c(names(gene_lists), "custom"), 
            inline = TRUE, 
            selected = "All")#,
          # disabled((selectInput("sel_custom_gene_set", label = "Custom gene lists", choices = "")))
        )
      })
      
      output$ui_sel_color_by  = renderUI({
        dt = meta_data()
        req(dt)
        ns <- session$ns
        n_unique = sapply(colnames(dt), function(nam){
          length(unique(dt[[nam]]))
        })
        facet_vars = names(n_unique)[n_unique <= 8]
        radioButtons(inputId = ns("sel_color_by"), label = "Color By", choices = facet_vars)
      })
      
      
      
      ## watch gene inputs
      observeEvent({
        input$sel_gene_list
        input$sel_custom_gene_set
        # input$txt_genes
      }, {
        sel = input$sel_gene_list
        if(sel == "custom"){
          gl = custom_gene_sets()[[input$sel_custom_gene_set]]
        }else{
          gl = NULL
          gl = gene_lists[[sel]]
          # message(paste(gl, collapse = ", "))
          gl = sort(unique(gl))
        }
        input_genes(gl)
      })
      
      output$ui_sel_sample_type_filter =  renderUI({
        req(meta_data())
        ns <- session$ns
        active_sample_types = unique(meta_data()$sample_type)
        checkboxGroupInput(
          inputId = ns("sel_sample_type_filter"), 
          label = "Samples Included", 
          choices = active_sample_types,
          selected = active_sample_types)
      })
      

      
      
      ##
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
      
      ## compare input genes with available expression data
      observeEvent({
        expression_data()
        input_genes()
      }, {
        req(expression_data())
        gl = input_genes()
        
        missed = setdiff(gl, rownames(expression_data()))
        if(length(missed) > 0){
          showNotification(paste("genes not present in expression:", paste(missed, collapse = ", ")), type = "warning")
        }
        valid_genes(setdiff(gl, missed))
      })
      ## notify about genes loaded
      observeEvent({
        valid_genes()
      }, {
        showNotification(paste0(length(valid_genes()), " genes loaded from ", input$sel_gene_list, "."))
      })
      
      observeEvent({
        custom_gene_sets()
      }, {
        gls = custom_gene_sets()
        if(length(gls) > 0) enable("sel_custom_gene_set")
        curr_sel = input$sel_custom_gene_set
        if(curr_sel %in% names(gls)){
          updateSelectInput(session, "sel_custom_gene_set", choices = names(gls), selected = curr_sel)
        }else{
          updateSelectInput(session, "sel_custom_gene_set", choices = names(gls))
        }
        showNotification(paste(length(gls), " custom gene sets"))
      })
      
      #### app_module_expression_matrix ####
      #running tsne
      # server_expression_matrix(
      #   input, output, session,
      #                          et,
      #                          active_dataset,
      #                          dataset_downstream,
      #                          expression_data,
      #                          meta_data,
      #                          tsne_input,
      #                          tsne_res)
      observeEvent({
        active_dataset()
      }, {
        et = active_dataset()
        expression_data(et$norm_counts)
        meta_data(et$meta_data)
        
        #reset downstream
        for(rv in dataset_downstream){
          rv(NULL)
        }
      })
      
      observe({
        showNotification(paste0("expression: ", nrow(expression_data()), " rows x ", ncol(expression_data()), " columns loaded."))
      })
      #### ####
      #### app_module_tsne ####
      # server_tsne(
      #   input, 
      #   output,
      #   session, 
      #   tsne_res, 
      #   tsne_input, 
      #   valid_genes, 
      #   meta_data, 
      #   code2type, 
      #   FACET_VAR)
      
      observeEvent({
        # expression_data()
        tsne_input()
        valid_genes()
        tsne_res()
      }, {
        showNotification("server_tsne")
        expr_mat = tsne_input()
        gl = valid_genes()
        req(expr_mat)
        req(gl)
        # req(is.null(tsne_res()))
        #choose dimensional reduction method, if possible
        if(ncol(expr_mat) < 3){
          showNotification("too_small")
          showNotification("Too few samples for dimensional reduction.", type = "error")
          tsne_worked = FALSE
        }else if(ncol(expr_mat) < 20){
          browser()
          showNotification("run_PCA")
          pc = prcomp(expr_mat[gl,])
          tsne_dt = as.data.table(pc$rotation[,1:2], keep.rownames = TRUE)[, c(2:3, 1)]
          setnames(tsne_dt, c("rn", "PC1", "PC2"), c("column_id", "tx", "ty"))
          tsne_worked = TRUE
        }else if(length(gl) > 0){
          tsne_worked = tryCatch({
            showNotification(paste("run_tsne:",  nrow(expr_mat[gl,]), "rows x", ncol(expr_mat[gl,]), "columns"))
            tsne_dt = run_TSNE(expr_mat[gl,])    
            TRUE
          }, error = function(e){
            FALSE
          })
        }else{
          tsne_worked = FALSE
        }
        
        if(tsne_worked){
          meta_dt = meta_data()
          tsne_dt = merge(tsne_dt, meta_dt, by = "column_id")
          tsne_dt[, tx := scales::rescale(tx, c(-.5, .5))]
          tsne_dt[, ty := scales::rescale(ty, c(-.5, .5))]
          tsne_res(tsne_dt) 
        }else{
          showNotification("Need more valid genes/samples to run t-sne.", type = "error")
          tsne_res(NULL)
        }
      })
      
      ### Plot t-sne
      output$plot_tsne <- renderPlot({
        req(tsne_res())
        req(input$sel_facet_var)
        browser()
        tsne_dt = tsne_res()
        tsne_dt = merge(tsne_dt, meta_data(), by = "column_id")
        if(input$sel_color_by ==  FACET_VAR$NONE){
          p = ggplot(tsne_dt, aes_string(x = "tx", y = "ty")) 
        }else{
          p = ggplot(tsne_dt, aes_string(x = "tx", y = "ty", color = input$sel_color_by))     
        }
        
        if(input$sel_facet_var != FACET_VAR$NONE){
          p = p + annotate("point", x= tsne_dt$tx, y = tsne_dt$ty, color = 'gray70', size = .3)
        }
        p = p + 
          geom_point() + 
          coord_fixed() +
          labs(x = "", y = "", title = "t-sne of TCGA samples", subtitle = paste(input$sel_gene_list, "gene list")) +
          theme(panel.background = element_blank(),
                panel.grid = element_blank())
        if(input$sel_facet_var != FACET_VAR$NONE){
          p = p + facet_wrap(paste0("~", input$sel_facet_var))
        }
        # }else if(input$sel_facet_var == FACET_VAR$SAMPLE_TYPE){
        #     p + facet_wrap(~sample_type)
        # }else if(input$sel_facet_var == FACET_VAR$PAM50){
        #     p + facet_wrap(~pam_call)
        # }else{
        #     stop("unrecognized input$sel_facet_var: ", input$sel_facet_var)
        # }
        p
      })
      #### ####
      
      #the gene expression mapped to tsne space
      server_gene_xy(input, output, session, tsne_res, expression_data, vis_gene, norm_description = et$norm_description)
      #upload user data for gene lists
      # gene_table, expression_data, custom_gene_sets
      server_upload(
        id = "default_upload",
        expression_data = expression_data, 
        custom_gene_sets = custom_gene_sets)
      #interface to select A and B set of points from scatterplot
      sample_groups = server_point_selection(id = "default_point_selection",
                                             tsne_clust = tsne_clust, 
                                             meta_data = meta_data, 
                                             tsne_res = tsne_res)
      sample_groups_A =sample_groups$A
      sample_groups_B =sample_groups$B
      
      #DESeq2
      DE_res = reactiveVal()
      
      observeEvent({
        sample_groups
      },{
        req(sample_groups)
        showNotification(paste0("A ", length(sample_groups_A()), "\n",
                                "B ", length(sample_groups_B())))
      })
      
      observeEvent({
        input$btn_runDE
      }, {
        req(tsne_input())
        req(sample_groups)
        diff_res = run_DE(tsne_input(), sample_groups_A(), sample_groups_B())
        DE_res(diff_res)
      })
      
      observeEvent({
        DE_res()
      }, {
        req(DE_res())
        showNotification(paste("diff gene count:", nrow(DE_res())))
      })
      
      output$dt_DE_res = DT::renderDataTable({
        if(is.null(DE_res())){
          DT::datatable(data.frame(waiting = "", for_ = "", data = ""))    
        }else{
          DT::datatable(DE_res())    
        }
        
      })
      
      #DEfast
      DE_fast_raw = reactiveVal()
      DE_fast_res = reactiveVal()
      
      observeEvent({
        input$btn_runDEfast
      }, {
        req(tsne_input())
        req(sample_groups)
        if(length(sample_groups_A()) == 0 | length(sample_groups_B()) == 0){
          DE_fast_raw(NULL)
          DE_fast_res(NULL)
        }else{
          dt = run_group.fast(tsne_input(), sample_groups_A(), sample_groups_B())
          DE_fast_raw(dt)
          p_dt = run_DE.fast(dt)   
          DE_fast_res(p_dt)    
        }
        
      })
      
      output$plot_boxes = renderPlot({
        req( DE_fast_raw())
        req(DE_fast_res())
        
        dt = DE_fast_raw()
        p_dt = DE_fast_res()
        
        p_dt
        high_dt = p_dt[lg2_min > 10][order(-abs(lg2_fc))][1:10][order(lg2_fc)]
        high_genes = as.character(high_dt$gene_name)
        
        sel_dt = dt[gene_name %in% high_genes]
        sel_dt$gene_name = factor(sel_dt$gene_name, levels = high_genes)
        ggplot(sel_dt, aes(x = gene_name, y = log2(expression+1), color = group)) +
          geom_boxplot(position = "dodge", width = .6) +
          labs(y = "log2 expression", x= "") +
          scale_color_manual(values = c("A" = "red", "B" = "blue", "." = "gray"), drop = FALSE) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      })
      
      output$plot_volcano = renderPlot({
        req(DE_fast_res())
        p_dt = DE_fast_res()
        
        p_dt
        high_dt = p_dt[lg2_min > 10][order(-abs(lg2_fc))][1:10]
        
        ggplot(p_dt, aes(x = lg2_fc, y = lg2_min)) +
          geom_point(alpha = .1) +
          geom_point(data= high_dt) +
          ggrepel::geom_text_repel(data = high_dt, aes(label = gene_name)) +
          labs(x = "log2 fold-change(B / A)", y = "log2 min(A, B)")
      })
    }
  )
}

#' expTSNE.runApp
#'
#' @param et expTSNE object
#'
#' @return
#' @export
#' @import shiny shinycssloaders
#' @rawNamespace import(shinyjs, except = runExample)
#' @examples
#' ex_data = system.file("extdata/test_expTSNE", package = "expTSNE", mustWork = TRUE)
#' et = expTSNE.load(ex_data)
#' expTSNE.runApp(et)
expTSNE.runApp = function(et, ...){
  # Define UI for application that draws a histogram
  ui <- fluidPage(
    theme = "bootstrap.css",
    # Application title
    useShinyjs(),
    titlePanel("TCGA t-sne"),
    expTSNE.ui_module("TSNE_panel")
  )
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
    expTSNE.server_module("TSNE_panel", et)
  }
  
  shiny::runApp(list(ui = ui, server = server), ...)
}
