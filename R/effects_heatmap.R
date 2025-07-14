#' @title Add Title Here
#'
#' @description This is a generic plotting utility for comparing the
#'   \dQuote{effects} on different features (rows) within a given
#'   signal (column), and for comparing the effect of the same
#'   feature across different signals.
#'
#' @param effects_matrix n x d numeric matrix, where n is the number
#'   of features and d is the number of different signals.
#'
#' @param zero_value Numbers smaller than \code{zero_value} (in
#'   magnitude) are not shown in the plot.
#'
#' @param size_range Passed as the \dQuote{range} argument to
#'   \code{\link[ggplot2]{scale_size}}.
#' 
#' @param font_size Specifies the font size for the plot.
#' 
#' @return Describe the return value here.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' set.seed(1)
#'X <- matrix(sample(seq(-5,5),200,replace = TRUE),20,10)
#' rownames(X) <- paste0("row",1:20)
#' colnames(X) <- paste0("col",1:10)
#' effects_heatmap(X,size_range = c(2,8),font_size = 9)
#' effects_heatmap(abs(X),size_range = c(2,8),font_size = 9)
#' 
#' @importFrom stats quantile
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_text
#' @importFrom cowplot theme_cowplot
#'
#' @export
#' 
effects_heatmap <- function (effects_matrix, zero_value = 0.01,
                             size_range = c(2,6), font_size = 10) {

  # Verify and process the effects_matrix input.
  if (!(is.matrix(effects_matrix) & is.numeric(effects_matrix)))
    stop("Input \"effects_matrix\" should be a numeric matrix")
  if (nrow(effects_matrix) < 2 | ncol(effects_matrix) < 2)
    stop("Input \"effects_matix\" should have at least 2 rows and ",
         "at least 2 columns")
  if (is.null(rownames(effects_matrix)) | is.null(colnames(effects_matrix)))
    stop("Input \"effects_matrix\" should be a named matrix; ",
         "that is, the rows and columns must be named. ",
         "See help(rownames) for details.")

  features <- rownames(effects_matrix)
  pdat <- data.frame(feature_name = features,stringsAsFactors = FALSE)
  pdat <- cbind(pdat,effects_matrix)
  pdat <- melt(pdat,id.vars = "feature_name",variable.name = "dim",
               value.name = "value")
  pdat <- transform(pdat,
                    effect_size = abs(value),
                    effect_sign = factor(sign(value),c(-1,0,1)),
                    feature_name = factor(feature_name,rev(features)))
  pdat$effect_size[abs(pdat$effect_size) < zero_value] <- NA
  effect_size_breaks <-
    unname(quantile(pdat$effect_size,probs = c(0,0.25,0.5,0.75,1),
                    na.rm = TRUE))
  if (any(pdat$effect_sign == -1))
    dot_colors <- c("navy","lightgray","orangered")
  else
    dot_colors <- c("slategray","lightgray","navy")
  return(ggplot(pdat,aes(x = dim,y = feature_name,size = effect_size,
                         fill = effect_sign)) +
         geom_point(color = "white",shape = 21,na.rm = TRUE) +
         scale_size(range = size_range,breaks = effect_size_breaks,
                    labels = round(effect_size_breaks,digits = 3)) +
         scale_fill_manual(values = dot_colors,drop = FALSE) +
         guides(size = guide_legend(override.aes = list(fill = "black")),
                fill = guide_legend(override.aes = list(color = "white",
                                                        shape = 21,
                                                        size = 3))) +
         labs(x = "",y = "",fill = "effect sign",size = "effect size") +
         theme_cowplot(font_size = font_size) +
         theme(panel.grid  = element_line(linewidth = 0.33,color = "gray"),
               axis.ticks  = element_blank(),
               axis.line   = element_blank(),
               axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
}
