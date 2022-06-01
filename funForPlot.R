g_legend<-function(a.gplot,direction = "horizontal"){
  a.gplot <- a.gplot +theme(legend.direction=direction)
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]] 
  return(legend)}


plotScore <- function(monocle,
                      traj,
                      sig,
                      shortName,
                      se = T,
                      ptBranching1 = NULL,
                      ptBranching2 = NULL) {
  
  
  plot <- ggplot(pData(monocle)[which(pData(monocle)$State %in% traj),] ,aes(color = AGE,y = get(sig), x=Pseudotime)) + 
    geom_smooth(se = se)  + 
    ylab("Score")+
    scale_color_manual(values = hue_pal()(2)[c(1,2)]) +
    ggtitle(shortName) + 
    theme(legend.text  = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 20))+
    guides(colour = guide_legend(override.aes = list(size=8)))
  
  if (!is.null(ptBranching1)) {
    plot <- plot + geom_vline(aes(xintercept=ptBranching1), color="red", linetype="dashed", size=1)
  }
  
  if (!is.null(ptBranching2)) {
    plot <- plot + geom_vline(aes(xintercept=ptBranching2), color="black", linetype="dashed", size=1)
  }
  
  return(plot)
}
