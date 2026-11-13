### check_colnames ###
#
#' Find required columns that are missing from a table.
#'
#' Returns the names in \code{cols} that are not columns of \code{x}. It does not stop on its own:
#' the calling function decides what to do, typically
#' \code{if (length(missing_cols) > 0) stop(...)}.
#'
#' @param x data.frame or sf object.
#' @param cols Character vector of required column names.
#'
#' @return Character vector of the missing column names; \code{character(0)} when none are missing.
#' @export
#'
#' @examples
#' check_colnames(data.frame(CATCHNUM = 1, network = "PB_0001"), c("CATCHNUM", "geometry"))
check_colnames <- function(x, cols) {
  setdiff(cols, colnames(x))
}
