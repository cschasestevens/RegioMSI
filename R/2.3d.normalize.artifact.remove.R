#### Remove artifact images from dataset ####

## Apply noise/artifact filter
fun.art.filter <- function(md,p.md) {
  if(class(md.samp[["ID"]]) != "factor") {
    md.samp[["ID"]] <- as.factor(md.samp[["ID"]])
    }
  
  d2 <- dplyr::left_join(
    md.samp,
    data.frame(
      d[c(1:p.md)],
      d[,names(d) %in% md.feat[["Input.name"]]]
    ),
    by = c("pixel","ID","Group")
  )
  
  ### Sort colnames and assign final names
  d2 <- data.frame(
    d2[,c(1:md)],
    d2[,sort(names(d2[,(md + 1):ncol(d2)]))]
  )
  
  d2.names <- dplyr::left_join(
    data.frame("Input.name" = names(d2[,(md + 1):ncol(d2)])),
    md.feat,
    by = "Input.name"
  )
  
  names(d2) <- c(names(d2[,c(1:md)]),d2.names[["Name.adduct"]])
  
  return(d2)
  
}




