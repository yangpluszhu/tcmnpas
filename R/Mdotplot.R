Mdotplot = function(object, x=~Cluster, colorBy="p.adjust", showCategory=5, by="geneRatio", includeAll=TRUE, font.size=12, title=""){ 
  df <- fortify__compareClusterResult(object, showCategory=showCategory, by=by, includeAll=includeAll) 
  Mplotting(df, x=x, type="dot", colorBy=colorBy, by=by, title=title, font.size=font.size) 
} 
 
