# Function to get the same colours as used by default in Seurat DimPlot
get_seurat_colours <- function(obj, group.by) {
	require(scales)
	identities <- obj@meta.data[,group.by]
	if (class(identities) != "factor") {
		identities <- factor(identities)
	}
	identities <- levels(identities)

	my_color_palette <- hue_pal()(length(identities))
	return(my_color_palette)
#usage:
# TSNEPlot(object = object, do.return = T) + 
# scale_color_manual(values = my_color_palette)

}


## Map Subtypes to general cell-type labels ##
simplify_annotations <- function(annotations, types=c("B", "Mac", "T", "Hep")) {
	simplified <- as.character(annotations)
	if ("B" %in% types) {
	simplified[simplified %in% c(
		"AntibodysecretingBcells", 
		"MatureBcells", "B_cells", "any_B_cell", "Anti_B")] <- "Bcells"
	}
	if ("Mac" %in% types) {
	simplified[simplified %in% c(
		"InflamatoryMacrophages", "Non-inflammatoryMacrophages")] <- "Macrophages"
	}
	if ("T" %in% types) {
	simplified[simplified %in% c(
		"CD3abTcells", "gdTcells1", "gdTcells2", "gd_T_1", "gd_T_2", "CD3_T")] <- "Tcells"
	}
	if ("Hep" %in% types) {
	simplified[simplified %in% c(
		"PericentralHep", "UnidentifiedHep", "PeriportalHep",
		"interzonalHep", "Hep", "PortalHep2", "PortaHep1", "PortHep", "CentralHep1", "CentralHep")] <- "Hepatocyte"
	}
	return(simplified);
}

Cell_type_colours <- rbind(
	c("ambiguous", "#7F7F7F"), #grey50
	c("doublet", "#7F7F7F"), #grey50
	c("Eryth", "#8B0000"), #darkred
	c("Hepatocyte", "#e41a1c"), #bright red
	c("OtherHep", "#e41a1c"), #bright red
	c("PortalHep", "#08306b"),
	c("interHep", "#08306b"),
	c("CentralHep", "#6baed6"),
	c("Stellate", "#fec44f"),
	c("Cholangiocyte", "#d94801"),
	c("Portalendo", "#993494"),
	c("Macrophage", "#377eb8"), # bright blue
	c("NonInfMac", "#2171b5"),
	c("InfMac", "#41b6c4"),
	c("Tcell", "#4daf4a"), #green
	c("gdTcells1", "#41ae76"),
	c("gdTcells2", "#00441b"),
	c("CD3abTcells", "#4daf4a"),
	c("Bcell", "#f781bf"), #violet
        c("AntiBcell", "#f78abf"), #violet
        c("MatBcell", "#6a51a3"), #violet
	c("NKcells", "#a65628"),
	c("LSECs", "#ffff33"),
	c("cvLSECs", "#fd8d3c"),
	c("PortalLSECs", "#4d004b")
)

colnames(Cell_type_colours) <- c("type", "colour")


map_cell_types <- function(types) {
	type <- as.character(types)
	tab <- unique(matrix(c("ambiguous", "ambiguous",
		"None", "ambiguous",
		"Unknown", "ambiguous",
		"Unknown1", "ambiguous",
		"Unknown2", "ambiguous",
		"Unknown3", "ambiguous",
		"Doublet", "doublet",
		"doublet", "doublet",
		"AntibodysecretingBcells", "AntiBcell",
		"AntiBcell", "AntiBcell",
		"Bcells", "Bcell",
		"B_cells", "Bcell",
		"any_B_cell", "Bcell",
		"Bcell", "Bcell",
		"CD3abTcells", "CD3abTcells",
		"CD3abTcell", "CD3abTcells",
		"NaiveTcell", "CD3abTcells",
		"CentralvenousLSECs", "cvLSECs",
		"CV_LSECs", "cvLSECs",
		"cvLSECs", "cvLSECs",
		"cvLSEC", "cvLSECs",
		"LSECs", "LSECs",
		"Cholangiocytes", "Cholangiocyte",
		"Cholangiocyte", "Cholangiocyte",
		"Cholan", "Cholangiocyte",
		"Chol", "Cholangiocyte",
		"Erythoidcells", "Eryth",
		"Erythoid", "Eryth",
		"Erythroid", "Eryth",
		"Erythroidcells", "Eryth",
		"Eryth", "Eryth",
		"gdTcells1", "gdTcells1",
		"gd_T_1", "gdTcells1",
		"gd_T_2", "gdTcells2",
		"gdTcell", "gdTcells2",
		"Hepatocyte", "Hepatocyte",
		"Hepatocytes", "Hepatocyte",
		"ErythHep", "Hepatocyte",
		"Hep", "Hepatocyte",
		"inflamatoryMacrophages", "InfMac",
		"InflamatoryMacrophages", "InfMac",
		"Monocyte-like", "InfMac",
		"InflammatoryMacrophages", "InfMac",
		"inflammatoryMacrophages", "InfMac",
		"InfMac", "InfMac",
		"Infl_Mac", "InfMac",
		"interzonalHep", "interHep",
		"interHep", "interHep",
		"InterHep", "interHep",
		"Macrophage", "Macrophage",
		"Mac", "Macrophage",
		"Macrophages", "Macrophage",
		"MatureBcells", "MatBcell",
		"MatBcell", "MatBcell",
		"NK-likecells", "NKcells",
		"NK_like", "NKcells",
		"NKTcell", "NKcells",
		"NKTcells", "NKcells",
		"NKcells", "NKcells",
		"Non-inflammatoryMacrophages", "NonInfMac",
		"NonInfl_Mac", "NonInfMac",
		"NonInflMac1", "NonInfMac",
		"NonInflMac22", "NonInfMac",
		"NonInfMac", "NonInfMac",
		"Kupffer", "NonInfMac",
		"PericentralHep", "CentralHep",
		"CentralHep", "CentralHep",
		"CentralHep1", "CentralHep",
		"PeriportalHep", "PortalHep",
		"PortaHep1", "PortalHep",
		"PortaHep2", "PortalHep",
		"PortalHep2", "PortalHep",
		"PortalHep1", "PortalHep",
		"PortalHep", "PortalHep",
		"PortHep", "PortalHep",
		"PeriportalLSECs", "PortalLSECs",
		"PeriLSECs", "PortalLSECs",
		"PortalLSECs", "PortalLSECs",
		"pLSECs", "PortalLSECs",
		"pLSEC", "PortalLSECs",
		"PortLSECs", "PortalLSECs",
		"Periportalendothelialcells", "Portalendo",
		"Portal_Endo", "Portalendo",
		"PortalEndo", "Portalendo",
		"Portal_Endo", "Portalendo",
		"PortEndo", "Portalendo",
		"Portalendothelialcells", "Portalendo",
		"Portalendo", "Portalendo",
		"Stellatecells", "Stellate",
		"Stellate", "Stellate",
		"Tcell", "Tcell",
		"Tcells", "Tcell",
		"gdTcells2", "gdTcells2",
		"gdTcells", "gdTcells2",
		"UnidentifiedHep", "OtherHep",
		"OtherHep", "OtherHep",
		"Bcell", "Bcell",
		"Cholan", "Cholangiocyte",
		"PortHep", "PortalHep",
		"Anti_B", "AntiBcell",
		"any_B_cell", "Bcell",
		"any_T_cell", "Tcell",
		"B_cells", "Bcell",
		"CD3_T", "CD3abTcells",
		"Bcells", "Bcell"), ncol=2, byrow=TRUE))

	colnames(tab) <- c("orig", "tidy");
	return(tab[match(types,tab[,1]),2])
}
	


type_2_colour <- function(type) {
	type <- as.character(map_cell_types(type))
	colour <- Cell_type_colours[match(type, Cell_type_colours[,"type"]),"colour"]
	return(colour)
}


Type_DimPlot <- function(myseur, type_col="consistent_labs", reduction="umap", cluster_col="seurat_clusters") {
	require(ggplot2)
	agg_coord_by_cluster <- function(coords, clusters) {
		x <- split(seq(nrow(coords)), clusters)
		xes <- sapply(x, function(a){median(coords[a,1])})
		yes <- sapply(x, function(a){median(coords[a,2])})
		return(rbind(xes, yes))
	}

	#umap_lab_pos <- agg_coord_by_cluster(Reductions(myseur, reduction)@cell.embeddings, myseur@meta.data[,cluster_col])
	umap_lab_pos <- agg_coord_by_cluster(myseur@reductions[[reduction]]@cell.embeddings, myseur@meta.data[,cluster_col])

	# UMAP + Ref scmap anno
	new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
	myseur@meta.data[,type_col] <- map_cell_types(myseur@meta.data[,type_col])
	new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% myseur@meta.data[,type_col],]

	DimPlot(myseur, reduction=reduction, group.by=type_col, pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=colnames(umap_lab_pos), colour="grey35")
}

convert_using_list <- function(vector_to_convert, conversion_list) {
       coverted <- sapply(vector_to_convert, function(x){conversion_list[[x]]})
       return(coverted)
}
