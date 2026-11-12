# gen_network_names
#' Create a vector of network names.
#'
#' Takes a vector of conservation area names and combines them into network names using the separator: \code{"__"}.
#'
#' Conservation area names should never include the separator: \code{"__"}, however \code{"_"} is acceptable.
#' 
#' By default all combinations of conservation area names will be created based on the \code{k} parameter. 
#' The length of the output will therefore equal \code{choose(length(in_names), k)}.
#' 
#' Input names are ordered in the network name based using \code{sort()}.
#' 
#' The output vector of names can be filtered using force_in to include specific names.
#'
#' @param in_names Vector of input names. Usually the names of conservation areas from the beaconsbuilder package (e.g. PB_0001) and/or user provided reserves (e.g. PA_1).
#' @param k Integer > 1 and <= \code{length(in_names)}. Count of names to appear in each output network name.
#' @param force_in Vector of names to force in to the output vector (usually a subset of in_names but could be a substring of an output 
#' name, see examples). Only network names containing one of these vectors will be included in the output.
#'
#' @return A vector of network names.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 2)
#' 
#' gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 2, "PB_1")
#' 
#' gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 3, c("PB_1", "PB_11"))
#' 
#' gen_network_names(c("PB_1", "PB_2", "PB_3", "PB_11"), 3, "PB_1__PB_11")

gen_network_names <- function(in_names, k, force_in = c()){
  
  # check k
  if(k < 2){
    stop("k should be 2 or more")
  }
  if(k > length(in_names)){
    stop("k cannot be > length(in_names)")
  }
  
  names_combined <- utils::combn(sort(as.character(in_names)), k, simplify=FALSE) # simplify=FALSE returns a list
  out_names <- sapply(names_combined, function(x) paste0(x,collapse="__"))
  
  if(length(force_in) > 0){
    out_names <- unique(unlist(lapply(force_in, function(x){
      c(
        out_names[grepl(paste0(x, "__"), out_names)], # matches the name followed by __
        out_names[grepl(paste0(x, "$"), out_names)] # matches the name when it ends the string
        # this regexp assumes no in_name is a complete duplicate of another in_name.
        # For example this will fail if there are two reserves named "a" and "aa" and you try to
        # force_in "a".
      )
    })))
  }
  return(out_names)
}

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