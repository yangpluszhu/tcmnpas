replaceDatabase = function(x){ 
  require(stringr) 
  temp=str_replace_all(x,'c\\(|\\)|"','') 
  return(temp) 
} 
 
