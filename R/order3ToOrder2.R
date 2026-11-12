order3ToOrder2 <- function(order3){
  # replicating builder function to get order 2 from 3
  base64Codes <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "=", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "_", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
  idx0 <- match(substr(order3, 1, 1), base64Codes) - 1 # C# indexes start at zero so minus 1
  idx1 <- match(substr(order3, 2, 2), base64Codes) - 1
  return(idx0 * length(base64Codes) + idx1)
}
