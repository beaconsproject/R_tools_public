# sep_network_names
#' Split network names into their constituent parts.
#'
#' Takes a vector of network names in the form output by [gen_network_names()] (i.e. \code{"__"} separated),
#' and splits the network names back into their original names.
#' 
#' If a single name is provided, a vector of names is returned. If a vector with length >1 is provided, a list
#' of vectors is returned with each list element named with the network name.
#'
#' @param network_names String of a network name, or a vector of multiple network names. Must include the \code{"__"} separator.
#'
#' @return A vector of names, or a list of vectors.
#' @export
#'
#' @examples
#' sep_network_names("PB_1__PB_2")
#' sep_network_names(c("PB_1__PB_2", "PB_1__PB_3"))
#'
#' x <- gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 2)
#' sep_network_names(x)
#'
#' x <- gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 4)
#' sep_network_names(x)

sep_network_names <- function(network_names){
  
  if(length(network_names) == 1){
    out_val <- strsplit(network_names, "__")[[1]]
  }
  
  if(length(network_names) > 1){
    out_val <- lapply(network_names, function(x){
      strsplit(x, "__")[[1]]
    })
    names(out_val) <- network_names
  }
  return(out_val)
}
