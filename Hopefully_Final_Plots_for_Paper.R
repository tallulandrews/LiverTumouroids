source("0_ColourScheme.R")



# Read in each manual clustering
expr_type <- "norm_exprs";
expr_mats <- list()
cell_colors <- list();
nn_graph <- list();
cell_coords <- list();

dimred_name <- "dm"

sce_objs <- list(CCA1="CCA1_manual_SC3.rds")

dim_reduction <- list(CCA1="CCA1_1000_Visualizations_dims.rds")

line_specific_genes <- list(CCA1=c("ATP1B3", "DPAGT1", "CLDN2", 
	"AQP5", "EZH2", "RECQL4", "TRAIP", "LMNB1", "CDCA7L", "IQGAP3", 
	"DDIAS", "CA9", "NDRG1", "ATP2B4", "DAPK1", "HIST1H2AC")) # genes for dotplots

heatmap_genes <- list(CCA1=c("CALM1", "DEGS2", "FASN", "FUT2", "MAP1LC3B", "ROIK3", 
	"HERPUD1", "EIF5AL1", "EIF5A", "CCT3", "HSPE1", "GOT2", "C1QBP", "LDHB", "MAD2L1",
	"ZWINT", "ASF1B", "CDK1", "RRM2", "NCAPH", "FEN1", "TYMS", "ANLN", "HMGB2", "SCD",
	"SCD", "NDRG1","ERO1A", "NDUFA4L2", "P4HA1", "QSOX1", "BNIP3L", "FXYD3"))

line_specific_groups <- list(CCA1=c("Progenitor", "Differentiated1", 
			"TICs", "Differentiated2")) # cluster names

scmap_results <- list(CCA1="CCA1_scmap_output.rds")


for (i in names(sce_objs)) {
	this_name <- paste(i, dimred_name, "Bespoke", sep="_")
	# save cluster IDs & markers
	# add lineage annotations to markers
	# do tSNE & sil_nn graph
	# save coords
	require("SingleCellExperiment")
	require("scater")
	require("Rtsne")
	require("M3Drop")
	require("CellTypeProfiles")
	require("Seurat")
	# Set up
	SCE <- readRDS(sce_objs[[i]])
	SCE <- SCE[!is.na(rowData(SCE)$biotype),]
	SCE <- SCE[rowData(SCE)$feature_symbol != "",]
	SCE <- SCE[!duplicated(rowData(SCE)$feature_symbol),]
	rownames(SCE) <- rowData(SCE)$feature_symbol 
	SCE <- SCE[rowData(SCE)$biotype == "protein_coding",]
	palette <- cluster_col(max(SCE$Manual_Clusters))
	cell_colours <- palette[SCE$Manual_Clusters]
	nCs <- factor_counts(SCE$Manual_Clusters)

	# save colours
	names(palette) <- 1:length(palette)
	SCE@metadata$palette <- palette
	SCE@metadata$C_keep <- nCs > 10
	SCE@metadata$C_names <- line_specific_groups[[i]]
	name_map <- SCE@metadata$C_names
	names(name_map) <- names(SCE@metadata$C_keep)[SCE@metadata$C_keep]


	coords <- readRDS(dim_reduction[[i]])
	
	coords <- coords[[dimred_name]]

	cell_keep <- SCE$Manual_Clusters %in% names(SCE@metadata$C_keep)[SCE@metadata$C_keep]
	SCE <- SCE[,cell_keep]
	SCE$named_clusters <- name_map[match(SCE$Manual_Clusters, names(name_map))]
	
	# ScatterPlot
	pdf(paste(this_name, "DimRedScatter.pdf", sep="_"), width=6, height=6)
	plot(coords$x[cell_keep], coords$y[cell_keep], col=cell_colours, 
		pch=16, xlab="DM_1", ylab="DM_2")
	dev.off()
	blank_plot <- function() {
		tmp <- par("mar")
		par(mar=c(0,0,0,0))
		plot(1,1, col="white", xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", main="", xlab="", ylab="", bty="n")
		par(mar=tmp)
	}
	pdf(paste(this_name, "DimRedScatter_Legend.pdf", sep="_"), width=6, height=6)
	blank_plot();
	my_order <- order(name_map, decreasing=T)
	legend("center", name_map[my_order], col=palette[my_order], pch=16, cex=2, bty="n")
	dev.off()


	# Marker - DotPlot
	require("ggplot2")
	seurat <- as.Seurat(SCE, data=expr_type)
	thing <- cbind(coords$x, coords$y)
	thing <- apply(thing, 2, function(x){
		x=x-min(x); x<-x/max(x)*2; x<-x-1})
	rownames(thing) <- colnames(SCE)
	colnames(thing) <- c("DM1", "DM2")
	seurat[["dm"]] <- CreateDimReducObject(embeddings = thing, key = "DM", assay = DefaultAssay(seurat))
	

	pdf(paste(this_name, "MarkerDot.pdf", sep="_"), width=8, height=5)
	a<-DotPlot(seurat, 
		features=line_specific_genes[[i]], 
		group.by="named_clusters")+
		theme(axis.text.x = element_text(angle = 90, hjust = 1))
	plot(a+labs(x="Genes", y="Type"))
	dev.off()
	# Marker - Scatters
	pdf(paste(this_name, "MarkerScatter.pdf", sep="_"), width=16, height=14)
	FeaturePlot(seurat, reduction="dm", features = line_specific_genes[[i]])
	dev.off()

	get_top_markers <- function(cluster, nmarks=5) {
		tmp <- FindMarkers(seurat, ident.1=cluster, group.by="named_clusters")
		tmp$detect_diff <- tmp$pct.1-tmp$pct.2
		tmp <- tmp[tmp$avg_logFC > 0 & tmp$detect_diff > 0,]
		return(rownames(tmp)[1:nmarks])
	}
	heatmap_genes <- c();
	for (n in unique(seurat@meta.data$named_clusters)) {
		heatmap_genes <- c(heatmap_genes, get_top_markers(n))
	}

	pdf(paste(this_name, "MarkerHeatmap.pdf", sep="_"), width=14, height=8)
	seurat <- ScaleData(seurat)
	my_order <- order(name_map, decreasing=F)
	DoHeatmap(seurat, features = c(heatmap_genes,line_specific_genes[[i]]),
		 size = 3, group.by = "named_clusters", group.colors=palette[my_order])
	dev.off()
	
	dat<- seurat@assays$RNA@data[rownames(seurat) %in% line_specific_genes[[i]],]
	dat <- data.frame(t(dat), seurat@meta.data$named_clusters)
	colnames(dat)[ncol(dat)] <- "type"
	
	# Marker Violin
	ggplot_palette <- palette
	names(ggplot_palette) <- name_map
	marker_violins <- list();
	for (gene in line_specific_genes[[i]]) {
		
	a <- ggplot(dat, aes_string(x="type", y=gene, fill="type"))+geom_violin()+scale_fill_manual(values=ggplot_palette)+ ggtitle(gene)+theme(plot.title = element_text(face="bold", size=30,hjust = 0.5 ), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))
	marker_violins[[gene]] <- a;
	}
	require(ggpubr)
	theme_set(theme_pubr())
	pdf(paste(this_name, "MarkerViolins.pdf", sep="_"), width=16, height=16)
	ggarrange(plotlist=marker_violins)
	dev.off()

	# ScatterPlot + scmap results
	scmap_out <- readRDS(scmap_results[[i]])
	pdf(paste(this_name, "ScmapScatter.pdf", sep="_"), width=6, height=6)
	plot(coords$x, coords$y, col="black", bg=scmap_out$scmap_cell_Cols, pch=21, xlab="DM_1", ylab="DM_2");
	lgend <- unique(scmap_out$scmap_cluster_labs);
	lgend_col <- unique(scmap_out$scmap_cell_Cols);
	dev.off()
}
