#' Formats DE from single cell into a matrix for heatmap
#'
#' @param folder  folder with text files with significant DE genes for single cell 
#' @param DElist list object with significant DE genes for single 
#' @param filetype file type in folder with significant DE genes
#' @param filecolsep separator for columns in files with significant DE genes
#' @keywords format single-cell data
#' @export
#' @import stringr
#' @examples
#' data <- format_data(folder='path/to/significant/DEresults/')
#' 

format_data <- function(folder=NULL, DElist=NULL, filetype=".txt", filecolsep="\t") {
    #This function formats data from DE csv files from LW output into a format that will then be passed '''
    #into generate figure function below that will generate heatmaps and dot plots of all significant data'''
    if (!is.null(folder)) {
        sig.files <-  list.files(path = folder, filetype, full.names = TRUE)
        sig.filesred <-  list.files(path = folder, filetype, full.names = FALSE)
        all.genes <- list()
        counter=1
        for (file in sig.files) {
            comparison <- str_remove(sig.filesred[counter], filetype)
            df <- read.table(file, sep=filecolsep, header=TRUE)
            df$genes <- rownames(df)
            df <- df[, c("avg_log2FC", "pct.1", "genes")]
            colnames(df)  <- c(paste0("avg_log2FC.", comparison), paste0("pct.1.", comparison), "genes")
            all.genes[[comparison]] <- df
            counter=counter+1
        }
        merged.all.genes <- Reduce(function(...) merge(..., by="genes", all = T), all.genes)
    } else if (!is.null(DElist)) {
        merged.all.genes <- Reduce(function(...) merge(..., by="genes", all = T), DElist)       
    } else if ((!is.null(folder) && !is.null(DElist))) {
        stop("Need either Folder or DElist notboth")
    }
    if (!is.null(folder)) {
        return(list("data"=merged.all.genes, "folders"=sig.files))
    } else {
        return(list("data"=merged.all.genes, "folders"=DElist))

    }
}
