# function for modifying asv tables with only taxa above certain threshold

## info input

# function to change taxon names if below set threshold
sort_abundant_taxa <- function(table, abundance_threshold, Abundance) {

  # check if rank names are correct
  if (identical(tail(colnames(phylo_melt), n=6), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))) {
    list_taxa_above_threshold <- function(input_table, rank, threshold, abundance_column) {
      df1 <- input_table %>% 
        group_by(across({{rank}})) %>% 
        summarise(max_abundance = max({{abundance_column}})) %>% 
        filter(max_abundance > threshold/100 & is.na({{rank}}) == FALSE) %>% 
        arrange(desc(max_abundance))
      out_vector <- pull(df1, {{rank}})
      return(out_vector)
    }
    
    ## create vectors for each rank which are above threshold
    genus_thr <- list_taxa_above_threshold(input_table = table,
                                           rank = Genus,
                                           threshold = abundance_threshold,
                                           abundance_column = Abundance)
    
    family_thr <- list_taxa_above_threshold(input_table = table,
                                            rank = Family,
                                            threshold = abundance_threshold,
                                            abundance_column = Abundance)
    
    order_thr <- list_taxa_above_threshold(input_table = table,
                                           rank = Order,
                                           threshold = abundance_threshold,
                                           abundance_column = Abundance)
    
    class_thr <- list_taxa_above_threshold(input_table = table,
                                           rank = Class,
                                           threshold = abundance_threshold,
                                           abundance_column = Abundance)
    
    phylum_thr <- list_taxa_above_threshold(input_table = table,
                                            rank = Phylum,
                                            threshold = abundance_threshold,
                                            abundance_column = Abundance)
    
    ## rename taxa, if taxa is on that rank in all samples below threshold, it will be named after other of the next higher rank
    pf1 <- table %>% 
      mutate(Genus_mod = if_else(Genus %in% genus_thr, Genus,
                             if_else(Family %in% family_thr, paste0("other_f_", Family, "_<", abundance_threshold, "%"), 
                                     if_else(Order %in% order_thr, paste0("other_o_", Order, "_<", abundance_threshold, "%"), 
                                             if_else(Class %in% class_thr, paste0("other_c_", Class, "_<", abundance_threshold, "%"), 
                                                     if_else(Phylum %in% phylum_thr, paste0("other_p_", Phylum, "_<", abundance_threshold, "%"), 
                                                             paste0("other_", Kingdom, "_<", abundance_threshold, "%")))))),
         Family_mod = if_else(Family %in% family_thr, Family, 
                              if_else(Order %in% order_thr, paste0("other_o_", Order, "_<", abundance_threshold, "%"), 
                                      if_else(Class %in% class_thr, paste0("other_c_", Class, "_<", abundance_threshold, "%"), 
                                              if_else(Phylum %in% phylum_thr, paste0("other_p_", Phylum, "_<", abundance_threshold, "%"), 
                                                      paste0("other_", Kingdom, "_<", abundance_threshold, "%"))))),
         Order_mod = if_else(Order %in% order_thr, Order, 
                             if_else(Class %in% class_thr, paste0("other_c_", Class, "_<", abundance_threshold, "%"), 
                                     if_else(Phylum %in% phylum_thr, paste0("other_p_", Phylum, "_<", abundance_threshold, "%"), 
                                             paste0("other_", Kingdom, "_<", abundance_threshold, "%")))),
         Class_mod = if_else(Class %in% class_thr, Class, 
                             if_else(Phylum %in% phylum_thr, paste0("other_p_", Phylum, "_<", abundance_threshold, "%"), 
                                     paste0("other_", Kingdom, "_<", abundance_threshold, "%"))),
         Phylum_mod = if_else(Phylum %in% phylum_thr, Phylum, 
                              paste0("other_", Kingdom, "_<", abundance_threshold, "%"))
         ) %>% 
      arrange(Phylum_mod, Class_mod, Order_mod, Family_mod, Genus_mod)
    
    # create list which contains individual vectors for each taxon rank
    output <- vector("list", 6)
    output[[1]] <- genus_thr
    output[[2]] <- family_thr
    output[[3]] <- order_thr
    output[[4]] <- class_thr
    output[[5]] <- phylum_thr
    # actual long output table
    output[[6]] <- pf1
    
    # give the elements of the list names which contains information of the threshold used
    names(output) <- c(paste0("genera_abv_", abundance_threshold, "_perc"),
                   paste0("families_abv_", abundance_threshold, "_perc"),
                   paste0("orders_abv_", abundance_threshold, "_perc"),
                   paste0("classes_abv_", abundance_threshold, "_perc"),
                   paste0("phyla_abv_", abundance_threshold, "_perc"),
                   paste0("ASV_table_taxa_abv_", abundance_threshold, "_perc"))
    
    # output list
    return(output)
  }
  else {print("Error: rank names have to be changed to 'Kingdom', 'Phylum', 'Class', 'Order', 'Family',  'Genus'")}
  
}