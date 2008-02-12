`A.chc` <-
function(hc) {
  n <- length(hc$order)
  hts <- diff(c(0,hc$hei))
  unA <- .C("UnA",as.integer(hc$mer),as.double(hts),as.integer(n),
          as.double(diag(n)),as.double(hts[1]*diag(n)),PACKAGE="compHclust")[[5]]
  return(matrix(unA,ncol=n)/sum(hts))
}

