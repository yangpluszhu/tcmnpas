tempcatch = function(x){ 
  require(stringr) 
  if (class(try(x$a$.attrs['href'],silent=T))=='try-error'){ 
    temp=x$a['href'] 
  }else{ 
    temp=x$a$.attrs['href'] 
  } 
  return(str_extract(temp,pattern='[0-9]+')) 
} 
 
