#' generates single cell heatmap 
#'
#' @param data data$data object from format_data function, to be run before this
#' @param colorder order that you want comparisons to show up on the heatmap
#' @param filename filename for heatmap (has to be .png, code will automatically generate pdf and svg versions of the heatmap)
#' @param cluster whether or not cluster genes in the heatmap
#' @param clusterbar whether to show clusters with a bar at the bottom of the heatmap
#' @param kmeanscluster whether to use kmeans clustering
#' @param kclusters Number of kmeans clusters if kmeanscluster=TRUE
#' @param distmethod distance method to use if kmeanscluster=FALSE (default=euclidean)
#' @param clustermethod cluster method to use default is set to ward.D2 but can uses any hclust method ("ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC))
#' @param cuttreeparam where to cut tree to generate clusters
#' @param figheight figure size height (inches)
#' @param figwidth figure size width (inches)
#' @param genefontsize size of gene names if clusterbar=FALSE
#' @param legendfont size of legend font
#' @param labelfont size of label (y axis) font
#' @param labels vector of specific labels
#' @param maxlfc maximum saturation level of LFC to show on heatmap
#' @param minlfc minimum saturation level of LFC to show on heatmap
#' @param linewidth maximum line width of individual genes markers in heatmap
#' @keywords format single-cell data
#' @export
#' @import plotrix
#' @import Polychrome
#' @import pheatmap
#' @import reshape2
#' @import textshape
#' @import ggplot2
#' @examples
#' clusters <- sc_heatmap(data$data)


sc_heatmap <- function(data, colorder=NULL, genes2remove=NULL,
    filename="heatmap.png", cluster=TRUE, clusterbar=TRUE, kmeanscluster=TRUE, kclusters=6, 
    distmethod="euclidean", clustermethod="ward.D2", cuttreeparam=6, figheight=3,figwidth=8, genefontsize=10, legendfont=6, labelfont=6,
    labels=NULL, maxlfc=NULL, minlfc=NULL,  linewidth=4) {
    set.seed(5)
    if (!is.null(genes2remove)) {
        data <- data[!(data$genes %in% genes2remove), ]
    }
    lfc <- c()
    pct <- c()
    for (i in colnames(data)) {
        if (grepl("avg_log2FC", i) == TRUE) {
            lfc <- c(lfc, i)
        } else if (i != "genes") {
            pct <- c(pct, i)
        }
    }
    datalfc <- data[, lfc]
    rownames(datalfc) <- data$genes
    datapct <- data[, pct]
    rownames(datapct) <- data$genes
    datalfc <- as.matrix(as.data.frame(datalfc))
    datalfc[is.na(datalfc)] <- 0
    datapct <- as.matrix(as.data.frame(datapct))
    class(datapct) <- "numeric"
    datapct[is.na(datapct)] <- 0

    colnames(datalfc) <- str_replace_all(colnames(datalfc), "avg_log2FC.", "")
    colnames(datapct) <- str_replace_all(colnames(datapct), "pct.1.", "")
    if (is.null(colorder)) {
        colorder=colnames(datalfc)
    }
    # some times column number to include will be different so make sure colorder is decided before clustering
    datalfc <- datalfc[, colorder]
    datapct <- datapct[rownames(datalfc), colorder]
    A <- rowSums(abs(datalfc))
    isexpr <- A > 0
    datalfc <- datalfc[isexpr, ]
    datapct <- datapct[isexpr, ]
    message("STATUS: number of DE genes ", length(rownames(datalfc)))
    if (isTRUE(cluster)) {
        if (isTRUE(kmeanscluster)) {
            hc <- kmeans(cbind(datalfc, datapct), kclusters)
            x <- sort(hc$cluster)
            clustergenes <- names(sort(hc$cluster))
        } else {
            if (distmethod != "pearson") {
                ds <- dist(cbind(datalfc, datapct), method = distmethod)
            } else {
                 ds <- as.dist(1 - cor(t(cbind(datalfc, datapct)), method = distmethod))
            }
            hc <- hclust(ds, method = clustermethod)
            x <- cutree(hc, h =cuttreeparam)
            geneorder <- hc$labels[hc$order]
            x <- x[geneorder]
            clustergenes <- names(x)
        }
        datalfccluster <- datalfc[clustergenes, colorder]
        dataclusterpct <- datapct[clustergenes, colorder]
    } else {
        datalfccluster <- datalfc[, colorder]
        dataclusterpct <- datapct[, colorder]
    }
    #values for color scheme need to be changed before I convert data to NAs 
    myPalette <- colorRampPalette(c("blue", "skyblue", "white", "orange", "red"))(100)
    if (!is.null(maxlfc)) {
        datalfccluster[datalfccluster > maxlfc] <- maxlfc
    }
    if (!is.null(minlfc)) {
        datalfccluster[datalfccluster < minlfc] <- minlfc
    }
    maxvalue <- max(datalfccluster)
    minvalue <- min(datalfccluster)
    sc6 <- scale_color_gradientn(
        colours = myPalette,
        values = scales::rescale(c(
            minvalue, 0,
            0, maxvalue
        ))
    )

    datalfccluster[datalfccluster == 0] <- NA
    dataclusterpct[dataclusterpct == 0] <- NA
    datalfcclustermelt <- reshape2::melt(datalfccluster)
    dataclusterpctmelt <- reshape2::melt(dataclusterpct)
    datalfcclustermelt$PCT.1 <- dataclusterpctmelt$value
    datalfcclustermelt$PCT.1 <- as.numeric(datalfcclustermelt$PCT.1)
    datalfcclustermelt$value <- as.numeric(datalfcclustermelt$value)
    colnames(datalfccluster) <- paste0("LFC.", colnames(datalfccluster))
    colnames(dataclusterpct) <- paste0("pct.", colnames(dataclusterpct))
    if (isTRUE(cluster)) {
        P36 <- createPalette(length(unique(x)), c("#ff0000", "#00ff00", "#0000ff"))
        names(P36) <- unique(x)
        colors <- c()
        colornames <- c()
        for (i in x) {
            colors <- c(colors, P36[[i]])
            v <- color.id(P36[[i]])
            colornames <- c(colornames, v[1])
        }
        names(colornames) <- clustergenes
        total <- merge(as.data.frame(colornames), datalfccluster, by = "row.names")
        colnames(total)[1] <- "genes"
        colnames(total)[2] <- "cluster"
        total <- merge(total, dataclusterpct, by.x = "genes", by.y = "row.names")
        rownames(total) <- total$genes
        total$genes <- NULL
        total <- total[rownames(datalfccluster), ] # order in same way as in heatmap
        write.csv(total, str_replace(filename, ".png", "totalinfo.csv"))
        write.csv(colornames, str_replace(filename, ".png", ".csv"))
        counter <- 0
        verticallines <- c()
        for (c in unique(colornames)) {
            verticallines <- c(verticallines, (counter + length(which(colornames == c))))
            counter <- (counter + length(which(colornames == c)))
        }
    } else {
        total <- merge(datalfccluster, dataclusterpct, by = "row.names")
        write.csv(total, str_replace(filename, ".png", "totalinfo.csv"))
    }

    datalfcclustermelt$Var2 <- factor(datalfcclustermelt$Var2, levels=colorder)
    pl <- ggplot(data = datalfcclustermelt, aes(
        x = Var1,
        y = Var2,
        color = value, size = PCT.1
    )) +
        geom_point(shape = "|") +
        scale_size_continuous(range = c(0, linewidth)) +
        sc6 + scale_y_discrete(labels = labels) +
        # scale_color_continuous_diverging(mid=0, l1=1.5) +
        labs(color = "LFCs", size = "% of cells exp.\ngene in Protected cond.") +
        theme_sc_heatmap() +
        theme(
            panel.grid.major.y = element_line(color = "white"),
            axis.text.y = element_text(size = labelfont),
            axis.title.y = element_blank(), axis.title.x = element_blank(),
            legend.title = element_text(size = legendfont), legend.text = element_text(size = legendfont), 
            legend.position = "top", legend.direction = "horizontal"
        ) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    if (isTRUE(cluster)) {
        if (isTRUE(clusterbar)) {
            pl + geom_vline(xintercept = verticallines, linetype = "dashed", color = "black", size = 0.1) + # nolint
            annotation_raster(raster = t(colors), 0, length(clustergenes), -1.0, 0.7, interpolate = T)
        } else { 
             pl + geom_vline(xintercept = verticallines, linetype = "dashed", color = "black", size = 0.1) + # nolint
            # scale_x_discrete(labels=) 
            theme(axis.text.x = element_text(size = genefontsize, angle = 90, hjust = 1, vjust = 0.5))
        }

    } else { 
        pl + theme(axis.text.x = element_text(size = genefontsize, angle = 90, hjust = 1, vjust = 0.5),
         axis.ticks.x = element_blank())
      }

    ggsave(filename, width = figwidth, height = figheight, dpi = 500, bg = "white")
    ggsave(str_replace(filename, ".png", ".svg"), width = figwidth, height = figheight, dpi = 500, bg = "white")
    ggsave(str_replace(filename, ".png", ".pdf"), width = figwidth, height = figheight, dpi = 500, bg = "white")
    if (isTRUE(cluster)) {
        return(colornames)
    } else {
        return(NULL)
    }
}


theme_sc_heatmap <- function(base_size = 14) {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.ticks.y = element_line(),
            axis.text = element_text(),
            axis.line.y = element_blank(), # element_line(colour = "black", size = .1),
            axis.line.x = element_blank(),
            axis.ticks = element_line(),
            panel.grid.major.y = element_line(colour = "black"),
            panel.grid.major.x = element_line(colour = "white"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(5, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}