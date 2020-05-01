
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
	tab <- matrix(c("ambiguous", "ambiguous",
		"AntibodysecretingBcells", "AntiBcell",
		"AntiBcell", "AntiBcell",
		"Bcells", "Bcell",
		"Bcell", "Bcell",
		"CD3abTcells", "CD3abTcells",
		"CentralvenousLSECs", "cvLSECs",
		"cvLSECs", "cvLSECs",
		"Cholangiocytes", "Cholangiocyte",
		"Cholangiocyte", "Cholangiocyte",
		"Erythoidcells", "Eryth",
		"Erythroidcells", "Eryth",
		"Eryth", "Eryth",
		"gdTcells1", "gdTcells1",
		"Hepatocyte", "Hepatocyte",
		"Hepatocytes", "Hepatocyte",
		"inflamatoryMacrophages", "InfMac",
		"InflamatoryMacrophages", "InfMac",
		"InflammatoryMacrophages", "InfMac",
		"inflammatoryMacrophages", "InfMac",
		"InfMac", "InfMac",
		"interzonalHep", "interHep",
		"interHep", "interHep",
		"Macrophage", "Macrophage",
		"Macrophages", "Macrophage",
		"MatureBcells", "MatBcell",
		"MatBcell", "MatBcell",
		"NK-likecells", "NKcells",
		"NKcells", "NKcells",
		"Non-inflammatoryMacrophages", "NonInfMac",
		"NonInfMac", "NonInfMac",
		"PericentralHep", "CentralHep",
		"CentralHep", "CentralHep",
		"PeriportalHep", "PortalHep",
		"PortalHep", "PortalHep",
		"PeriportalLSECs", "PortalLSECs",
		"PortalLSECs", "PortalLSECs",
		"Periportalendothelialcells", "Portalendo",
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
		"Bcells", "Bcell"), ncol=2, byrow=TRUE)

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
	
