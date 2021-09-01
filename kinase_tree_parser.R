# Kinase tree parser functions
library(ggplot2)
library(plyr)
library(dplyr)

mapData = function(mappings, data, szScale = list(from = c(0,2), to = c(0, 50)),
                   clScale = list(low = -1, mid = 0, high = 1),
                   clValues = c("#3bef3b", "#000000", "#f42411")) {
  if (class(mappings) == "formula") {
    fterms   <- terms(mappings)
    mappings <- rownames(attr(fterms, "factors"))
  }
  df       <- data.frame(names = data[[mappings[1]]], clrdata = data[[mappings[2]]], szdata = data[[mappings[3]]])
  df$size  <- rescale(df$szdata, from = szScale$from, to = szScale$to, clip = TRUE)
  df$color <- mapToDivColorMap(df$clrdata, clScale, clValues)
  df
}

kinhubParse = function(mappings, data, szScale =list(from = c(0,2), to = c(0, 50)),
                       clScale = list(low = -1, mid = 0, high = 1),
                       clValues = c("#3bef3b", "#000000", "#f42411") ){
  
  df <- mapData(mappings, data, szScale, clScale, clValues)
  parsed_str <- dlply(df, ~names, .fun = function(x) {
    crgb <- col2rgb(x$color)
    str  <- c(paste("@0:", x$size,":",crgb[1], ",", crgb[2], ",", crgb[3]), as.character(x$names))
  })
  return(parsed_str)
}

treeFile = function(fname, mappings, data, parsefun = kinhubParse, ...){
  pstr = parsefun(mappings, data, ...)
  fid = file(description = fname, open = "wt")
  lapply(pstr, FUN = writeLines, con = fid)
  close(fid)
}

rescale = function (x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE), clip = FALSE) {
  if (zero_range(from) || zero_range(to)) {
    return(ifelse(is.na(x), NA, mean(to)))
  }
  x <- (x - from[1])/diff(from) * diff(to) + to[1]
  
  if (clip) {
    # elegant clipping by Josh O'Brien @ Stackoverflow
    #http://stackoverflow.com/questions/13868963/clip-values-between-a-minimum-and-maximum-allowed-value-in-r
    x <- to[1] + (x > to[1])*(x-to[1])  - (x > to[2]) * (x - to[2])
  }
  x
}

mapToDivColorMap = function(x, clScale = list(low = -1, mid = 0, high = 1), clValues = c("#3bef3b", "#000000", "#f42411")) {
  scaled.color <- rescale(x - clScale$mid, from = c(clScale$low, clScale$high), clip = TRUE)
  aScale       <- scale_color_gradient2(low = clValues[1], mid = clValues[2], high = clValues[3])
  aScale$palette(scaled.color)
}