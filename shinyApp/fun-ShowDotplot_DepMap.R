ShowDotplot_DepMap <- function(DepMap.info.all, gMCS.info.all, gmcs_database, gene.target.info, by_gMCS = F,
                               database, database_filter_mode, database_filter_selected, threshold_value,
                               flag_database_filter_show_only, flag_LinearRegression = T, flag_show_SYMBOL = T,
                               flag_color_by_partner_gene = F, database_unit){
  
  # browser()
  
  # remove empty database
  # database_filter_selected <- database_filter_selected[database_filter_selected != "---"]
  
  # extract the neccesary info from the target
  gene.ENSEMBL <- gene.target.info$ENSEMBL
  gene.SYMBOL <- gene.target.info$SYMBOL
  
  
  # test the data
  if(gene.ENSEMBL %in% gMCS.info.all$gMCSs.ENSEMBL.txt) { return(tibble(x = 1, y = 1) %>% ggplot(aes(x, y)) + theme_void() + geom_text(label ="This is an essential gene, it has no partner")) }
  # if(sum(gMCS.info.all$gMCSs.ENSEMBL.max[,gene.ENSEMBL])==1){by_gMCS = F; flag_color_by_partner_gene = F}
  
  if (by_gMCS) {
    gmcs.ENSEMBL <- strsplit(gMCS.info.all$gMCSs.ENSEMBL.txt[as.numeric(as.character(gene.target.info$gMCS))],'--')[[1]]
    gmcs.ENSEMBL <- setdiff(gmcs.ENSEMBL, gene.ENSEMBL)
    gmcs.SYMBOL <- strsplit(gMCS.info.all$gMCSs.SYMBOL.txt[as.numeric(as.character(gene.target.info$gMCS))],'--')[[1]]
    gmcs.SYMBOL <- setdiff(gmcs.SYMBOL, gene.SYMBOL)
    
    sdf <- reshape2::dcast(DepMap.info.all$DepMapGeneExpression %>% filter(UNIT==database_unit), ENSEMBL_target~DepMap_ID, value.var = "logTPM", fun.aggregate = sum)
    sdf <- sdf %>% filter(ENSEMBL_target %in% gmcs.ENSEMBL)
    sdf <- reshape2::melt(apply(sdf[,-1],2,max))
    sdf <- data.frame(DepMap_ID = rownames(sdf),
                      logTPM = sdf$value,
                      ENSEMBL_target = gene.ENSEMBL,
                      UNIT = database_unit)
    
  } else {
    sdf <- DepMap.info.all$DepMapExpressionByGene %>% 
      filter(ENSEMBL_target==gene.ENSEMBL) %>% 
      filter(gMCS.database==gmcs_database) %>% 
      filter(UNIT==database_unit) 
  }
  
  # browser()
  
  if (flag_color_by_partner_gene & !by_gMCS) {
    
    
    if (flag_show_SYMBOL) { 
      sdf <- merge(sdf, 
                   gMCS.info.all$table.genes.HumanGEM %>% 
                     rename(ENSEMBL_partner = ENSEMBL) %>% 
                     rename(SYMBOL_partner = SYMBOL)  %>% 
                     rename(ENTREZID_partner = ENTREZID))
      sdf$partner_gene <- sdf$SYMBOL_partner
    } else {
      sdf$partner_gene <- sdf$ENSEMBL_partner
    }
    
    # reorder partner genes, set a maximum of four and add percentage 
    # (only show those with more than 1% of samples)
    browser
    zz <- sort(table(sdf$partner_gene),decreasing = T)
    zz_levels <- zz_labels <- names(zz)
    zz_per <- as.numeric(zz)/sum(zz)*100
    nn <- min(sum(zz_per>1),3)
    if (length(zz_labels)>nn){
      nn <- 1:(nn)
      zz_labels[-nn] <- "others"
      zz_per[-nn] <- sum(zz_per[-nn])
    }
    zz_labels <- paste0(zz_labels, " (", round(zz_per,1),"%)" )
    sdf$partner_gene <- factor(as.character(sdf$partner_gene),
                               levels = zz_levels, labels = zz_labels)
    
  } else {
    sdf$partner_gene <- "noColor"
  }
  
  # head(sdf)
  
  
  sdf.2 <- DepMap.info.all$DepMapEssentialityByGene %>% 
    filter(ENSEMBL_target==gene.ENSEMBL) %>% 
    filter(essentiality.database==database)
  
  sdf <- merge(sdf, sdf.2)
  sdf$logTPM <- as.numeric(as.character(sdf$logTPM))
  sdf$essentiality_score <- as.numeric(as.character(sdf$essentiality_score))
  
  # head(sdf)
  
  # add the filter selected
  sdf$isFilterSelected <- "rest"
  
  if (database_filter_mode!="none"){  # & database_filter_selected!="---"){
    sdf <- merge(sdf, DepMap.info.all$dictionary.CCLE[,c("DepMap_ID", database_filter_mode)])
    sdf$isFilterSelected[sdf[,database_filter_mode] %in% database_filter_selected] <- sdf[sdf[,database_filter_mode] %in% database_filter_selected, database_filter_mode]
  }
  sdf$isFilterSelected <- as.factor(sdf$isFilterSelected)
  
  if (database_unit == "zscores(TPM+1))") {
    database_unit = "zscores(TPM))"
  } else if (database_unit == "proteomics") {
    database_unit = "Relative Protein Expression"
  }

  
  XLAB <- paste0("Minimum expression of associated gMCSs partner genes \\[ ", sub("log2", "log_2", database_unit), " \\]")
  if (by_gMCS){
    if (flag_show_SYMBOL) { 
      XLAB <- paste0("Maximum expression of all partner genes: \\textit{", paste(gmcs.SYMBOL, collapse = ", "), "} \\[ ", sub("log2", "log_2", database_unit), " \\]")
    } else {
      XLAB <- paste0("Maximum expression of all partner genes: \\textit{", paste(gmcs.ENSEMBL, collapse = ", "), "} \\[ ", sub("log2", "log_2", database_unit), " \\]")
    }
  }
  
  
  # set the Y lab
  if (flag_show_SYMBOL) { 
    YLAB <- paste0("Essentiality of \\textit{",gene.SYMBOL , "} \\[ ",database," score \\]")
  } else {
    YLAB <- paste0("Essentiality of \\textit{",gene.ENSEMBL , "} \\[ ",database," score \\]")
  }
  
  ## Generate a Dotplot
  pp_dotplot <- ggplot(sdf, aes(x = logTPM, y = essentiality_score)) + 
    # scale_y_continuous() + 
    theme_classic() + 
    theme(text = element_text(size=12),
          axis.text.x = element_text(color='black', size=12),
          # legend.title = element_blank(),
          # axis.ticks.x = element_blank(),
          axis.text.y = element_text(color='black', size=12),
          legend.text = element_text(color='black', size=12,  margin = margin(r = 30, unit = "pt")),
          # legend.spacing.x = unit(3, 'cm'),
          # legend.margin = margin(t = 0, r = 2, b = 1, l = 2, unit = "cm"),
          plot.title = element_text(color='black', size=20),
          # axis.line.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          strip.text = element_text(color='black', size=14), 
          axis.title = element_text(color = "black", size = 14, face = "bold"),
          legend.position = ifelse(sdf$partner_gene[1]=="noColor", "none", "right")) + 
    xlab(latex2exp::TeX(XLAB)) + ylab(latex2exp::TeX(YLAB)) 
  
  if (!is.na(threshold_value)) {
    pp_dotplot <- pp_dotplot + geom_hline(yintercept = threshold_value, linetype = "dotdash", col = "red3", size = 1)
  }
  
  if (database_filter_mode!="none"){
    if (flag_database_filter_show_only){
      sdf3 <- sdf[sdf$isFilterSelected!="rest",]
      pp_dotplot <- pp_dotplot + 
        geom_point(mapping = aes(col = partner_gene, shape = isFilterSelected), size = 2, data = sdf3) +
        # scale_color_manual(values = c("grey50", "black"), breaks = c("0","1")) +
        scale_shape_manual(values = c(21, rep(16, length(database_filter_selected))), breaks = c("rest",database_filter_selected))  
      # scale_size_manual(values = c(1.5, 3), breaks = c("0","1"))
    } else {
      sdf3 <- sdf[sdf$isFilterSelected!="rest",]
      pp_dotplot <- pp_dotplot + 
        geom_point(mapping = aes(col = partner_gene, shape = isFilterSelected), size = 1.5, data = sdf) +
        geom_point(mapping = aes(col = partner_gene, shape = isFilterSelected), size = 3, data = sdf3) +
        # scale_color_manual(values = c("grey50", "black"), breaks = c("0","1")) +
        scale_shape_manual(values = c(21, rep(8, length(database_filter_selected))), breaks = c("rest",database_filter_selected))  
      # scale_size_manual(values = c(1.5, 3), breaks = c("0","1"))
    }
  } else {
    pp_dotplot <- pp_dotplot + geom_point(aes(col = partner_gene), size = 2)
  }
  
  if (flag_LinearRegression){
    if (flag_database_filter_show_only){
      sdf3 <- sdf[sdf$isFilterSelected!="rest",]
    } else { sdf3 <- sdf}
    pp_dotplot <- pp_dotplot +  geom_smooth(method = "lm", color = "black", fullrange = T,
                                            mapping = aes(x = logTPM, y = essentiality_score), formula = 'y ~ x',
                                            data = sdf3)
    
    z <- cor.test(sdf3$essentiality_score, sdf3$logTPM)
    pp_dotplot <- pp_dotplot + annotate("text", label = paste0("\u03C1 = ",round(z$estimate,3),"           \np-value = ",format.pval(z$p.value),"           "),
                                        x = Inf, y = -Inf, hjust = 1, vjust = -1, size = 4)
  }
  
  if (flag_color_by_partner_gene){
    pp_dotplot <- pp_dotplot + scale_color_brewer(type = "qual", palette = "Dark2") + 
      guides(color=guide_legend(title="Partner gene / biomarker"))
  } else {
    pp_dotplot <- pp_dotplot + scale_color_manual(values = c("black"), breaks = c("noColor"))
  }
  
  # pp_dotplot
  return(pp_dotplot)
}
