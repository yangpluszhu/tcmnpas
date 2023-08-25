theme_dose = function(font.size=14){ 
  theme_bw() %+replace% 
    theme(axis.text.x = element_text(colour = "black", 
                                     size = font.size, vjust =1 ), 
          axis.title.x = element_text(colour="black", 
                                      size = font.size), 
          axis.text.y = element_text(colour = "black", 
                                     size = font.size, hjust =1 ), 
          axis.title.y = element_text(colour="black", 
                                      size = font.size, angle=90) 
    ) 
} 
 
