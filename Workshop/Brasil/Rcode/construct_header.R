construct_header <- function(df, grp_names, span, align = "c", draw_line = T) {
  if (length(align) == 1) align <- rep(align, length(grp_names))
  if (!all.equal(length(grp_names), length(span), length(align)))
    stop("grp_names and span have to have the same length!")
  
  if (ncol(df) < sum(span)) stop("Span has to be less or equal to the number of columns of df") 
  
  header <- mapply(function(s, a, grp) sprintf("\\multicolumn{%i}{%s}{%s}", s, a, grp),
                   span, align, grp_names)
  header <- paste(header, collapse = " & ")
  header <- paste0(header, " \\\\")
  
  if (draw_line) {
    # where do we span the lines:
    min_vals <- c(1, 1 + cumsum(span)[1:(length(span) - 1)])
    max_vals <- cumsum(span)
    line <- ifelse(grp_names == "", "", 
                   sprintf("\\cmidrule(lr){%i-%i}", min_vals, max_vals))
    line <- paste(line[line != ""], collapse = " ")
    
    header <- paste0(header, "  ", line, "\n  ")
  }
  
  addtorow <- list(pos = list(-1, -1, nrow(df)),
                   command = c("\\hline\n  ", header, "\\hline\n  "))
  return(addtorow)
}