
server_expression_matrix = function(input, output, session,
                                    app_datasets,
                                    active_dataset,
                                    dataset_downstream,
                                    tcga_data,
                                    meta_data,
                                    tsne_input,
                                    tsne_res
){
    #handle change in dataset selection
    observeEvent({
        active_dataset()
    }, {
        et = active_dataset()
        tcga_data(et$norm_counts)
        meta_data(et$meta_data)
        
        #reset downstream
        for(rv in dataset_downstream){
            rv(NULL)
        }
    })
    
    observe({
        showNotification(paste0("expression: ", nrow(tcga_data()), " rows x ", ncol(tcga_data()), " columns loaded."))
    })
    ### Sample Type
    # observeEvent({
    #     input$sel_sample_type_filter
    #     tcga_data()
    # }, {
    #     req(tcga_data())
    #     samp_dt = sample_data()
    #     sel_types = input$sel_sample_type_filter
    # 
    #     samp_dt = samp_dt[sample_type %in% sel_types]
    #     k = colnames(tcga_data()) %in% samp_dt$sample_id
    #     tsne_input(tcga_data()[, k])
    #     tsne_res(NULL)
    # })
}