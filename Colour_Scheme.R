
Cell_type_colours <- rbind(
	c("ambiguous", "#7F7F7F"), #grey50
	c("Eryth", "#8B0000"), #darkred
	c("Hepatocyte", "#e41a1c"), #bright red
	c("OtherHep", "#e41a1c"), #bright red
	c("PortalHep", "#084081"),
	c("interHep", "#0868ac"),
	c("CentralHep", "#2b8cbe"),
	c("Stellate", "#fec44f"),
	c("Cholangiocytes", "#ec7014"),
	c("Portalendo", "#993494"),
	c("Macrophage", "#377eb8"), # bright blue
	c("NonInfMac", "#377eb8"),
	c("InfMac", "#41b6c4"),
	c("Tcell", "#4daf4a"), #green
	c("gdTcells1", "#41ae76"),
	c("gdTcells2", "#00441b"),
	c("CD3abTcells", "#4daf4a"),
	c("Bcell", "#f781bf"), #violet
	c("AntiBcell", "#756bb1"), #violet
	c("MatBcell", "#bcbddc"), #violet
	c("NKcells", "#a65628"),
	c("LSECs", "#ffff33"),
	c("cvLSECs", "#ff7f00"),
	c("PortalLSECs", "#4d004b"))


map_cell_types <- function(types) {
	type <- as.character(types)
	tab <- unique(matrix(c("ambiguous", "ambiguous",
		"AntibodysecretingBcells", "AntiBcell",
		"AntiBcell", "AntiBcell",
		"Bcells", "Bcell",
		"Bcell", "Bcell",
		"CD3abTcells", "CD3abTcells",
		"CentralvenousLSECs", "cvLSECs",
		"CV_LSECs", "cvLSECs",
		"cvLSECs", "cvLSECs",
		"LSECs", "LSECs",
		"Cholangiocytes", "Cholangiocyte",
		"Cholangiocyte", "Cholangiocyte",
		"Erythoidcells", "Eryth",
		"Erythoid", "Eryth",
		"Erythroidcells", "Eryth",
		"Eryth", "Eryth",
		"gdTcells1", "gdTcells1",
		"gd_T_1", "gdTcells1",
		"gd_T_2", "gdTcells2",
		"Hepatocyte", "Hepatocyte",
		"Hepatocytes", "Hepatocyte",
		"Hep", "Hepatocyte",
		"inflamatoryMacrophages", "InfMac",
		"InflamatoryMacrophages", "InfMac",
		"InflammatoryMacrophages", "InfMac",
		"inflammatoryMacrophages", "InfMac",
		"InfMac", "InfMac",
		"Infl_Mac", "InfMac",
		"interzonalHep", "interHep",
		"interHep", "interHep",
		"Macrophage", "Macrophage",
		"Mac", "Macrophage",
		"Macrophages", "Macrophage",
		"MatureBcells", "MatBcell",
		"MatBcell", "MatBcell",
		"NK-likecells", "NKcells",
		"NK_like", "NKcells",
		"NKcells", "NKcells",
		"Non-inflammatoryMacrophages", "NonInfMac",
		"NonInfl_Mac", "NonInfMac",
		"NonInfMac", "NonInfMac",
		"PericentralHep", "CentralHep",
		"CentralHep", "CentralHep",
		"CentralHep1", "CentralHep",
		"PeriportalHep", "PortalHep",
		"PortaHep1", "PortalHep",
		"PortaHep2", "PortalHep",
		"PortalHep", "PortalHep",
		"PortHep", "PortalHep",
		"PeriportalLSECs", "PortalLSECs",
		"PeriLSECs", "PortalLSECs",
		"PortalLSECs", "PortalLSECs",
		"Periportalendothelialcells", "Portalendo",
		"Portal_Endo", "Portalendo",
		"Portalendothelialcells", "Portalendo",
		"Portalendo", "Portalendo",
		"Stellatecells", "Stellate",
		"Stellate", "Stellate",
		"Tcell", "Tcell",
		"Tcells", "Tcell",
		"gdTcells2", "gdTcells2",
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
	


Cell_type_colours <- rbind(
	c("ambiguous", "#7F7F7F"), #grey50
	c("Eryth", "#8B0000"), #darkred
	c("Hepatocyte", "#e41a1c"), #bright red
	c("OtherHep", "#e41a1c"), #bright red
	c("PortalHep", "#084081"),
	c("interHep", "#0868ac"),
	c("CentralHep", "#2b8cbe"),
	c("Stellate", "#fec44f"),
	c("Cholangiocyte", "#ec7014"),
	c("Portalendo", "#993494"),
	c("Macrophage", "#377eb8"), # bright blue
	c("NonInfMac", "#377eb8"),
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
	c("cvLSECs", "#ff7f00"),
	c("PortalLSECs", "#4d004b"))
colnames(Cell_type_colours) <- c("type", "colour")



type_2_colour <- function(type) {
	type <- as.character(map_cell_types(type))
	colour <- Cell_type_colours[match(type, Cell_type_colours[,"type"]),"colour"]
	return(colour)
}
	
simplify_annotations <- function(annotations, types=c("B", "Mac", "T", "Hep")) {
	simplified <- as.character(annotations)
	if ("B" %in% types) {
	simplified[simplified %in% c(
		"AntibodysecretingBcells", 
		"MatureBcells")] <- "Bcells"
	}
	if ("Mac" %in% types) {
	simplified[simplified %in% c(
		"InflamatoryMacrophages", "Non-inflammatoryMacrophages")] <- "Macrophages"
	}
	if ("T" %in% types) {
	simplified[simplified %in% c(
		"CD3abTcells", "gdTcells1", "gdTcells2")] <- "Tcells"
	}
	if ("Hep" %in% types) {
	simplified[simplified %in% c(
		"PericentralHep", "UnidentifiedHep", "PeriportalHep",
		"interzonalHep")] <- "Hepatocyte"
	}
	return(simplified);
}
