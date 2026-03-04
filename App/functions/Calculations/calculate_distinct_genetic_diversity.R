calculate_DGD <- function(input, reactive_DGD_table, reactive_KBA_info) {
  req(input$file)
  file_ext <- tools::file_ext(input$file$name)
  
  # Laden der Daten basierend auf der Dateiendung
  if (file_ext == "csv") {
    df <- read.csv(input$file$datapath)
    pop <- as.vector(df[, 2])
    genind <- df2genind(df[2:ncol(df)], pop = pop, sep = ", ")
  } else if (file_ext %in% c("gen", "dat", "gtx")) {
    genind <- import2genind(input$file$datapath)
  } else if (file_ext %in% c("str", "stru")) {
    req(input$col_area, input$n_ind)
    genind <- import2genind(input$file$datapath,
                            onerowperind = FALSE, 
                            n.ind = input$n_ind, 
                            n.loc = input$n_loc, 
                            col.lab = input$col_geno, 
                            col.pop = input$col_area, 
                            ask = FALSE)
  } else {
    showNotification("Invalid file format. Please upload a valid file.", type = "error")
    return(NULL)
  }
  
  # Berechnungen durchführen
  # split genind: 
  genind_list <- lapply(unique(genind@pop), function(subset_by_potential_KBAs) {
    individuals <- which(genind@pop == subset_by_potential_KBAs)
    genind[individuals, ]
  })
  names(genind_list) <- unique(genind@pop)
  
  # create a function to calculate distinct genetic diversity
  calculate_EDplus <- function(genind) {
    allele_frequencies <- tab(genind, freq = TRUE)
    dist_matrix <- bitwise.dist(genind)
    mod <- taxondive(t(allele_frequencies), dist_matrix)
    return(mod$EDplus)
  }
  
  # apply the function to every genind in our list
  EDplus_values <- lapply(genind_list, calculate_EDplus)
  

  # Ergebnisse berechnen
  sumDplus <- sum(unlist(EDplus_values))
  ratios <- unlist(EDplus_values)/sumDplus * 100
  
  # Ergebnisse in reactive Datenrahmen speichern
  reactive_DGD_table$df_ratios <- data.frame(
    site = names(ratios),
    `distinct genetic diversity (Δ<sup>+</sup>)` = unlist(EDplus_values),
    `distinct genetic diversity (Δ<sup>+</sup>) [%]` = ratios,
    check.names = FALSE
  )
  
  # KBA-Informationen aktualisieren
  reactive_KBA_info$A1a <- ratios[ratios >= 0.5]
  reactive_KBA_info$A1b <- ratios[ratios >= 1]
  reactive_KBA_info$A1c <- ratios[ratios >= 0.1]
  reactive_KBA_info$A1d <- ratios[ratios >= 0.2]
  reactive_KBA_info$B1  <- ratios[ratios >= 10]
}