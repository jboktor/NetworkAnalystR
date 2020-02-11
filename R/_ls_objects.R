###########
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) BiocGenerics::sapply(names, function(x)
    fn(AnnotationDbi::get(x, pos = pos)))
  names <- AnnotationDbi::ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) rlang::as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- data.table::ifelse(rlang::is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    backports::capture.output(timeDate::format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- BiocGenerics::t(napply(names, function(x)
    h2o::as.numeric(dim(x))[1:2]))
  vec <- rlang::is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- dplyr::data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  ipred::print(BiocGenerics::lapply(dataSet, object.size));
  globals::names(out) <- S4Vectors::c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!gdata::missing(order.by))
    out <- out[BiocGenerics::order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- git2r::head(out, n)
  out
}
