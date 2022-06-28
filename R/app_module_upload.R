


ui_upload = function(id = "default_upload"){
  ns <- NS(id)
  ex_dir = system.file(package = "expTSNE", "extdata/example_data", mustWork = TRUE)
  ex_files = dir(ex_dir, full.names = TRUE)
  names(ex_files) = basename(ex_files)
  
  sidebarLayout(
    sidebarPanel(
      
      fluidRow(
        column(textInput(inputId = ns("txt_gene_set_name"), label = "Gene set name", value = "custom"), width = 9),
        column(actionButton(inputId = ns("btn_add_gene_set"), label = "Add"), width = 3, style = "margin-top: 24px;",)
      ),
      tags$hr(),
      tabsetPanel(
        id = ns("gene_list_method"),
        tabPanel("Paste", value = "paste",
                 textAreaInput(ns("txt_genes"), label = "Custom Genes", value = "paste genes here")#,
        ),
        tabPanel("Upload", value = "upload",
                 fileInput(inputId = ns("BtnUploadPeakfile"), label = "Browse Local Files"),
                 actionButton(inputId = ns("btnExampleData"), label = "Load Selected Data"),
                 selectInput(inputId = ns("selExampleData"), label = "Select Example Data", choices = ex_files)    
                 
        )
      ),
      uiOutput(outputId = ns("rpt_custom_gene_sets")),
      tags$h5("TODO: gene set import/export"),
      tags$h5("TODO: filtering after upload"),
      tags$h5("TODO: hide example data")
    ),
    mainPanel(
      DT::dataTableOutput(outputId = ns("DT_PasteGenes_DataFrame"))
    )
  )
  
}

server_upload = function(id = "default_upload", expression_data, custom_gene_sets = NULL){
  ex_dir = system.file(package = "expTSNE", "extdata/example_data", mustWork = TRUE)
  ex_files = dir(ex_dir, full.names = TRUE)
  names(ex_files) = basename(ex_files)
  
  if(is.null(custom_gene_sets)){
    custom_gene_sets = reactiveVal(list())
  }
  moduleServer(
    id,
    function(input, output, session){
      ### File Upload
      #Set Preview reactives
      PreviewSet_Filepath = reactiveVal(value = "", label = "PreviewSet_Filepath")
      PreviewSet_Name = reactiveVal(value = "", label = "PreviewSet_Name")
      PreviewSet_DataFrame = reactiveVal()
      gene_table = reactiveVal(data.frame())
      
      observeEvent({
        PreviewSet_Filepath()
        PreviewSet_Name()
      }, {
        if(PreviewSet_Filepath() == ""){
          PreviewSet_DataFrame(NULL)
        }else{
          # showNotification(paste0("loading ", PreviewSet_Name()))
          out = decide_parse_FUN(PreviewSet_Filepath(), PreviewSet_Name())
          PreviewSet_DataFrame(out)
        }
      })
      
      observeEvent(input$BtnUploadPeakfile, {
        PreviewSet_Filepath(input$BtnUploadPeakfile$datapath)
        PreviewSet_Name(input$BtnUploadPeakfile$name)
      })
      
      observeEvent(input$BtnCancelFile, {
        PreviewSet_Filepath("")
        PreviewSet_Name("")
      })
      
      observeEvent(input$btnExampleData,
                   {
                     PreviewSet_Filepath(input$selExampleData)
                     PreviewSet_Name(names(ex_files)[which(input$selExampleData == ex_files)])
                   })
      
      observeEvent(
        PreviewSet_DataFrame(),
        {
          req(PreviewSet_DataFrame())
        }
      )
      
      observe({
        if(input$gene_list_method == "upload"){
          # showNotification("DEBUG list by upload")
          df = PreviewSet_DataFrame()[[1]]
          if(is.null(df)){
            gene_table(data.frame())  
          }else if(nrow(df) == 0){
            gene_table(data.frame())
          }else{
            gene_table(df)
          }
        }
      })
      
      ### Pasting genes
      #genes parsed for paste/upload
      parsed_genes = reactiveVal()
      observeEvent({
        input$txt_genes
      }, {
        # showNotification(paste("parsed_genes update", id))
        gl = parse_gl(input$txt_genes)
        parsed_genes(gl)
      })
      
      observe({
        if(input$gene_list_method == "paste"){
          # showNotification(paste("gene_table update", id))
          gl = parsed_genes()
          if(is.null(gl)){
            gene_table(data.frame())  
          }else if(length(gl) == 0){
            gene_table(data.frame())  
          }else{
            gene_table(data.frame(gene_name = gl))
          }    
        }
        
      })
      display_table = reactiveVal()
      observeEvent({
        gene_table()
      }, {
        # showNotification(paste("display_table update", id))
        df = gene_table()
        df = as.data.frame(df)
        valid_genes = rownames(expression_data())
        if(nrow(df) > 0){
          df = locate_genes_in_df(df, valid_genes)
          df = validate_genes_in_df(df, valid_genes)    
        }else{
          df = data.table(A = "please upload data")
        }
        display_table(df)
      })
      
      ### Show table
      output$DT_PasteGenes_DataFrame = DT::renderDataTable({
        DT::datatable(display_table())    
      })
      
      ### Report valid
      
      ### Add gene set
      observeEvent({
        input$btn_add_gene_set
      }, {
        current = custom_gene_sets()
        df = display_table()
        gl = subset(df, valid == TRUE)$gene_name
        valid_genes = rownames(expression_data())
        gl = merge(data.frame(gl = gl, upper = toupper(gl)), 
                   data.frame(valid = valid_genes, upper = toupper(valid_genes)), 
                   by = 'upper')$valid 
        current[[input$txt_gene_set_name]] = gl
        custom_gene_sets(current)
      })
      
      output$rpt_custom_gene_sets = renderUI({
        ns = session$ns
        gls =  custom_gene_sets()
        hl = lapply(names(gls), function(nam){
          x = gls[[nam]]
          tags$tr(
            tags$td(nam, style = "padding:0 15px 0 3px;"),
            tags$td(length(x))
          )
        })
        tags$table(#cellpadding = "10", cellspacing = "10",
          tags$tr(
            tags$th("Name", style = "padding:0 15px 0 3px;"),
            tags$th("Length")
          ),
          hl
        )
      })
      
      return(list(reactive_custom_gene_sets = custom_gene_sets))
    })
}
