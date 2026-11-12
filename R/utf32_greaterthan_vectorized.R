utf32_greaterthan_vectorized <- function(x, y){
  sapply(x, function(xx){
    utf32_greaterthan(xx, y)
  })
}
