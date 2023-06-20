plot_2d <- function (PrjStruct, SeuratObj, reduction = 'dc',
                     dims = c(1,2), color_cells_by = "orig.ident", 
                     color_palette = NULL, viridis_option = 'D',
                     density_show = F, padding = 0.1, n_bins = 1000, bw = 0.2, density_cutoff = 0.3, 
                     cell_size = 2, alpha = 1, node_show = F,node_size = 1, node_label_show = F,
                     trajectory_color = "black", trajectory_width = 1, trajectory_alpha = .7)
{
  
  x <- dims[[1]]
  y <- dims[[2]]
  
  

  Embeddings <- Embeddings(SeuratObj,reduction = reduction)
  Metadata <- SeuratObj@meta.data %>% rownames_to_column(var = 'sample_name')
  Expression <- GetAssayData(SeuratObj, assay = DefaultAssay(SeuratObj), slot = 'data') %>% as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column(var = 'sample_name')
  
  data_df <- Embeddings[,c(x,y)] %>% as.data.frame()
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df <- data_df %>% rownames_to_column(var = 'sample_name')
  data_df <- data_df %>% left_join(Metadata, by = 'sample_name') %>% left_join(Expression, by = 'sample_name')
  
  xlims <- c(min(data_df[,'data_dim_1']),max(data_df[,'data_dim_1']))
  ylims <- c(min(data_df[,'data_dim_2']),max(data_df[,'data_dim_2']))
# cellの位置関係をplot------------------------------------------------------------------------
  plot_df <- data_df %>% dplyr::select(cell_color = color_cells_by, everything())
  
  if(methods::is(plot_df$cell_color, "numeric"))
  {
    p <- ggplot(data = plot_df) +
      geom_point(aes(x = data_dim_1, y = data_dim_2, color = cell_color), size = I(cell_size), alpha = I(alpha)) +
      ggraph::scale_color_viridis(name = color_cells_by, option = viridis_option) +
      xlab(paste0(toupper(reduction),'_',x)) + ylab(paste0(toupper(reduction),'_',y)) +
      theme_classic()
  }
  else {
    if(is.null(color_palette)){
      color_palette <- scales::hue_pal()(n = length(unique(plot_df$cell_color)))
      names(color_palette) <- unique(plot_df$cell_color)
    }
    p <- ggplot(data = plot_df) +
      geom_point(aes(x = data_dim_1, y = data_dim_2, color = cell_color), size = I(cell_size), alpha = I(alpha)) +
      scale_color_manual(name = color_cells_by, values = color_palette, breaks = names(color_palette)) + 
      xlab(paste0(toupper(reduction),'_',x)) + ylab(paste0(toupper(reduction),'_',y)) +
      theme_classic()
  }
  
# groupを囲む枠をつくる------------------------------------------------------------------------
  if(density_show){
  xlims <- c(min(data_df$data_dim_1), max(data_df$data_dim_1))
  ylims <- c(min(data_df$data_dim_2), max(data_df$data_dim_2))
  
  xpad <- diff(xlims) * padding
  ypad <- diff(ylims) * padding
  
  xlims <- xlims + c(-xpad, xpad)
  ylims <- ylims + c(-ypad, ypad)
  
  xbw <- diff(xlims) * bw
  ybw <- diff(ylims) * bw
  
 group_density <- data_df %>%
    dplyr::select(group_id = color_cells_by, data_dim_1, data_dim_2) %>%
    nest(positions = c(data_dim_1, data_dim_2)) %>% 
    mutate(contour = map2(positions, group_id, function(positions, group_id) {
      density <- MASS::kde2d(positions$data_dim_1, positions$data_dim_2, h = c(xbw, ybw), lims = c(xlims, ylims), n = n_bins)
      level <- max(density$z) * density_cutoff
      contour <- with(density, contourLines(x, y, z, levels = level))
      map2_df(contour, seq_along(contour), function(contour, contour_i) {
        tibble(
          data_dim_1 = contour$x,
          data_dim_2 = contour$y,
          contour_id = paste0(group_id, "_", contour_i)
        )
      })
    })) %>%
    unnest(contour)
  
  polygon <- geom_polygon(data = group_density, aes(x = data_dim_1, y = data_dim_2, fill = group_id, group = contour_id), alpha = 0.4)
  if(is.null(color_palette)){
    p <- p + polygon
  }else{
    p <- p + polygon + scale_fill_manual(name = color_cells_by, values = color_palette)
  }
  }
  
# Trajectoryのnodeの位置関係を抽出----------------------------------------------------------------
  ica_space_df <- PrjStruct$NodePositions %>% as.data.frame() %>% 
    dplyr::select(prin_graph_dim_1 = x, 
                  prin_graph_dim_2 = y) %>% 
    mutate(node_num = seq(1,nrow(.)))
  dp_mst <- PrjStruct$Edges %>% as.data.frame() %>% dplyr::select(from = 1, to = 2)
  edge_df <- dp_mst %>% 
    dplyr::left_join(ica_space_df %>% 
                       dplyr::select(from = "node_num", 
                                     from_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     from_prin_graph_dim_2 = "prin_graph_dim_2"), by = "from") %>% 
    dplyr::left_join(ica_space_df %>% 
                       dplyr::select(to = "node_num", 
                                     to_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     to_prin_graph_dim_2 = "prin_graph_dim_2"), by = "to") %>% 
    mutate(edge.id = row_number())
  
  if(node_show){
    p <- p +
      ggraph::geom_edge_link(data = edge_df, 
                             width = trajectory_width, color = trajectory_color, alpha = trajectory_alpha,
                             aes(x = from_prin_graph_dim_1, xend = to_prin_graph_dim_1, 
                                 y = from_prin_graph_dim_2, yend = to_prin_graph_dim_2)) +
      geom_point(data = ica_space_df, aes(x = prin_graph_dim_1, y = prin_graph_dim_2)) 
    
  }else{
    p <- p +
      ggraph::geom_edge_link(data = edge_df, 
                             width = trajectory_width, color = trajectory_color, alpha = trajectory_alpha,
                             aes(x = from_prin_graph_dim_1, xend = to_prin_graph_dim_1, 
                                 y = from_prin_graph_dim_2, yend = to_prin_graph_dim_2)) 
  }
  
  if(node_label_show){
    p <- p +
      ggrepel::geom_text_repel(data = ica_space_df, aes(label = node_num,x = prin_graph_dim_1, y = prin_graph_dim_2))
  }
  
  
  

  p + xlim(xlims) + ylim(ylims)
}
　