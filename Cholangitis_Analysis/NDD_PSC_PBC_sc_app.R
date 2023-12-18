#### Packages

library(shiny)
set.seed(101)

library(ggplot2)
library(Matrix)
library(scales)
require(dplyr)
library(cowplot)



# Read data
input_data <- readRDS("data/SC_Integrated_Map_absolute_min_shinydata.rds")

# Create a select list coloured based on significance
all_genes <- rownames(input_data$mean)
is.sig_PSC <- apply(input_data$FDR_PSC, 1, min) < 0.05
is.sig_PBC <- apply(input_data$FDR_PBC, 1, min) < 0.05
colour.sig <- "<span style='color:black';>"
colour.no.sig <- "<span style='color:red';>"
gene_select <- rep(colour.no.sig, length(is.sig_PSC))
gene_select[is.sig_PSC | is.sig_PBC] <- colour.sig
gene_select <- paste(gene_select, all_genes, "</span>", sep="")
names(all_genes) <- gene_select

# Separate cell-types into groups to make plot more readable
cell_types <- list(
		"Myeloid"=c("ActMac", "Kupffer--LSEC-Doublet", "Kupffer", "LAM-like", "MHCII", "Monocyte", "NKT--Mac-Doublet", "pDC", "cDC", "Neutrophil", "Hepato--Mac"),
		"NKT"=c("NKT--Mac-Doublet", "CD4T", "CD3T-lrNK", "CD8T-cNK", "CD8T", "cNK", "lrNK", "cvLSEC--T-Doublet", "MAST", "MatB--CD4T-Doublet", "MatB", "AntiB", "NKT", "Tcell", "Prolif"),
		"Parenchymal"=c("Fibroblast", "Arterial", "C-Hepato2", "C-Hepato", "CholMucus", "Chol", "cvEndo", "cvLSEC", "cvLSEC--T-Doublet", "Hepato", "I-Hepato", "Kupffer--LSEC-Doublet", "P-Hepato2", "P-Hepato", "ppLSEC", "Stellate", "RBC", "Hepato--Mac"))


# Get data to plot on UMAP
all_umap_features <- colnames(input_data$metadata)[c(1:3,6:7)]
umap_default = "cell_type"

# defaults
no_labels_symbol <- "(None)"
default_umap <- "cell_type"
default_genes <- c("MARCO", "CD5L", "ACP5", "TREM2", "LYZ", "VCAN", "IL1B")
default_genes_symbol <- "(Example GeneSet)"
current_genes <- c();
default_typeclass <- "Myeloid"

# Colours
get_colours <- function(identities) {
        require(scales)
        if (class(identities) != "factor") {
                identities <- factor(identities)
        }
        identities <- levels(identities)
        my_color_palette <- hue_pal()(length(identities))
	  names(my_color_palette) <- identities
        return(my_color_palette)
}

# UMAP Plot

agg_coord_by_cluster <- function(coords, clusters) {
      x <- split(seq(nrow(coords)), clusters)
      result <- sapply(x, function(a) {apply(coords[a,],2,median)})
      return(result)
}

Dim_plot <- function(input_data, feature, type_subset=NULL, label=TRUE) {
	this_palette <- get_colours(input_data$metadata[,feature])
	if (feature == "cell_type") {

	}
	thisplot <- ggplot(input_data$metadata, aes(UMAP1, UMAP2))+ geom_point(aes(color=.data[[feature]]), size=0.2) + theme_classic()

	if (feature %in% c("clusters_res_0.5", "cell_type")) {
		lab_pos <- agg_coord_by_cluster( as.matrix( input_data$metadata[,c("UMAP1", "UMAP2")] ), input_data$metadata[,feature] )
		group_labels <- colnames(lab_pos)
		if (sum(group_labels %in% type_subset) > 1) {
			# Only show cell-type labels for the cell-types in the DotPlot
			group_labels[!group_labels %in% type_subset] <- ""
		}
		thisplot <- thisplot + annotate("text", x=lab_pos[1,], y=lab_pos[2,], label=group_labels, color="black", fontface = "bold") + 
				guides(color="none")
	}
	return(thisplot)
}



# DotPlot
My_DotPlot_shiny <- function (my_mean_expr_mat, my_detection_rate_mat, features, type_class, # Input
				cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, # Colours of Dots
    				dot.min = 0, dot.scale = 6, scale = TRUE, scale.by = "radius", scale.min = 0.1, scale.max = 5) { # Size of Dots

	# Select Data
	columns <- colnames(my_mean_expr_mat)
	columns <- columns[!columns %in% c("CD4T--RBC-Doublet_NDD", "CD4T--RBC-Doublet_PSC", "cNK--RBC-Doublet_PSC", 
							"Hepato--Mac_NDD", "Hepato--Mac_PBC", "Hepato--Mac_PSC", 
							"Mac--B-Doublet_NDD", "Mac--B-Doublet_PBC", "Mac--B-Doublet_PSC",
							"Mac--Fibro-Doublet_NDD", "Mac--Fibro-Doublet_PBC", "Mac--Fibro-Doublet_PSC",
							"MatB--RBC_PBC", "MatB--RBC_PSC", "RBC_NDD", "RBC_PBC", "RBC_PSC", "NKT_PSC")]

	phenotypes_col <- unlist(lapply(strsplit(columns, "_"), function(x){x[[2]]}))
	types_col <- unlist(lapply(strsplit(columns, "_"), function(x){x[[1]]}))
	columns <- columns[types_col %in% cell_types[[type_class]]]

	features <- features[features %in% rownames(my_mean_expr_mat)]
	avg.exp <- my_mean_expr_mat[features,columns]/input_data$scale.factor
	pct.exp <- my_detection_rate_mat[features,columns]

	phenotypes_col <- unlist(lapply(strsplit(columns, "_"), function(x){x[[2]]}))
	types_col <- unlist(lapply(strsplit(columns, "_"), function(x){x[[1]]}))
	
	# Turn into dataframe for ggplot
	if (length(features)==1) {
		data.plot <- data.frame(avg.exp=as.vector(avg.exp), pct.exp=as.vector(pct.exp), 
					features.plot= features, 
					id=names(avg.exp),
					phenotype=phenotypes_col)
	} else {
		data.plot <- data.frame(avg.exp=as.vector(avg.exp), pct.exp=as.vector(pct.exp), 
					features.plot= rep(rownames(avg.exp), times=ncol(avg.exp)), 
					id=rep(colnames(avg.exp), each=nrow(avg.exp)),
					phenotype=rep(phenotypes_col, each=nrow(avg.exp)) )
	}
	#rownames(data.plot) <- data.plot$features.plot # DUPLICATE ROWNAMES
	#data.plot$id <- factor(data.plot$id, lavels=cell_types[[type_class]]) # THIS DOESN'T WORK

	scale.func <- switch(EXPR = scale.by, size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))

	# Scale expression for each gene
	avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == 
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
		    data.use[data.use < col.min] <- col.min
		    data.use[data.use > col.max] <- col.max
            }
            else {
                data.use <- log1p(x = data.use)
            }
            return(data.use)
        })
	
	avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
	data.plot$features.plot <- factor(x = data.plot$features.plot, 
        levels = features)
	data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
	data.plot$pct.exp <- data.plot$pct.exp



	# Turn expression into Colours 
	Null.colour="lightgrey"; 
	NDD.colour = colorRampPalette(colors=c(Null.colour, "red"))(20); 
	PBC.colour = colorRampPalette(colors=c(Null.colour, "purple"))(20); 
	PSC.colour = colorRampPalette(colors=c(Null.colour, "blue"))(20); 

	avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
	data.plot$avg.exp.scaled <- avg.exp.scaled
	dot.colour <- apply(data.plot, 1, FUN=function(x) {
		if (x[5]=="NDD") {return(NDD.colour[as.numeric(x[6])])}
		if (x[5]=="PSC") {return(PSC.colour[as.numeric(x[6])])}
		if (x[5]=="PBC") {return(PBC.colour[as.numeric(x[6])])}
		})
	data.plot$colors <- dot.colour
      color.by <- "colors"
	colour_scheme <- c(NDD.colour, PSC.colour, PBC.colour)
	names(colour_scheme) <- colour_scheme

	# Determing Dot Size
	data.plot[data.plot$pct.exp < 1/1000, "pct.exp"] <- 1/1000
#	if (!is.na(x = scale.min)) {
#		data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
#	}
#	if (!is.na(x = scale.max)) {
#		data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
#	}

	# Add significance
	data.plot$cell_type <- sapply(strsplit(data.plot$id, "_"), function(x){x[[1]]})
	pscFDR <- input_data$FDR_PSC[as.character(data.plot$features.plot), types_col]
	pscFDR <- apply(data.plot, 1, FUN=function(x) {pscFDR[x[3],x[8]] })
	pbcFDR <- input_data$FDR_PBC[as.character(data.plot$features.plot), types_col]
	pbcFDR <- apply(data.plot, 1, FUN=function(x) {pbcFDR[x[3],x[8]] })
	data.plot$pscFDR <- pscFDR
	data.plot$pbcFDR <- pbcFDR
	data.plot$FDR <- 1
	data.plot$FDR[data.plot$phenotype == "PSC"] <- data.plot$pscFDR[data.plot$phenotype == "PSC"]
	data.plot$FDR[data.plot$phenotype == "PBC"] <- data.plot$pbcFDR[data.plot$phenotype == "PBC"]
	data.plot$label = ""
	data.plot$label[data.plot$FDR < 0.05] <- "*"
#	data.plot$label[data.plot$FDR < 0.005] = "**"
#	data.plot$label[data.plot$FDR < 0.0005] = "***"

	# Make the plot
	plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
		y = "id", label="label")) + geom_point(mapping = aes_string(size = "pct.exp", 
		color = color.by)) + theme_cowplot() + scale_color_manual(values = colour_scheme, guide="none")+
		geom_text(hjust=0.5, vjust=+0.65, size=6, colour="green", fontface="bold") +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		theme(axis.title.x=element_blank(), axis.title.y=element_blank())
	return(plot)
}





#==============================================================



#### UI
ui <- fluidPage(
	titlePanel("Single Cell RNAseq of NDD, PSC, PBC livers"),

	# Sidebar Layout with input and output definitions ----
	sidebarLayout (

		# Sidebar panel for inputs ----
		sidebarPanel(

			# Input: Select Gene
		  
    			selectizeInput("gene", "Choose a gene (red = not significantly DE):", choices = all_genes, 
                  	options = list(render = I("{
					item: function(item, escape) { return '<div>' + item.label + '</div>'; },
					option: function(item, escape) { return '<div>' + item.label + '</div>'; }
				}"))), 

			# Buttons:
			# Add
			actionButton("addgene", "Add"),
			# Clear
			actionButton("cleargenes", "Clear"),

			# Input: Remove Gene
			#selectInput(inputId = "to_rm",
	            #      label = "Choose a gene to remove:",
                  #	choices = c(" ", current_genes), 
			#	selected = " "),
			uiOutput('to_rm'),

			# Remove
			actionButton("removegene", "Remove"),
			

			# Cluster Subset
			selectInput(inputId = "subset", 
				label = "Display Subset:",
				choices=names(cell_types),
				selected=c("Myeloid")),

			# UMAP feature
			selectInput(inputId = "umapfeature", 
				label = "feature for UMAP:",
				choices=all_umap_features,
				selected=c("cell_type")),

			checkboxInput(inputId="labelumap", "Label UMAP?", value = FALSE, width = NULL),

			### Downloading Plot ###
			# Input: Set resolution of image to download
			sliderInput("res",
				label = "Download DotPlot Resolution (dpi)",
                       	min = 50, 
				max = 300, 
				step = 50,
				value = 150),
	

			# Download Button
			downloadButton("downloadPlot", "Download DotPlot")
			
		),
	
		# Main panel for displaying outputs ----
		mainPanel(

			# Formatted Text for explanation & Citation information
			strong("Citation:"), 
			"Andrews, T.S.*, Atif, J.* Liu, J.C.*, Perciani, C.T., Ma, X.Z., Thoeni, C., Slyper, M., Eraslan, G., Segerstolpe, A., Manuel, J., Chung, S., Winter, E., Cirlan, I., Khuu, N., Fischer, S., Rozenblatt-Rosen, O., Regev, A., McGilvray, I.E., Bader, G.D., MacParland, S.A.,", 
			"(2021)", em("Hepatology Communications."), strong("doi :"), "10.1002/hep4.1854",
			tags$a(href="https://doi.org/10.1002/hep4.1854", "Link to Paper"),

			# Geneset
			textOutput(outputId = "Geneset"),

			# Output: UMAP
			plotOutput(outputId = "UMAPPlot"),
			strong("Figure 1 (above):"), "UMAP coloured by your feature of choice. If showing cell-types, only cell-types included in the Dot Plot below are labelled.",

			strong("Figure 2 (below):"), "Dot Plot of the gene expression of the selected genes across the selected subset of cell-types. 
			Red = Healthy (NDD, n=6), Purple = Primary Billiary Cholangitis (PBC, n=2), Blue = Primary Sclerosing Cholangitis (PSC, n=8)
			Green stars indicate a gene was significantly different (FDR < 5%) from the NDD control.",
			# Output : Dot Plot
			plotOutput(outputId = "DotPlot")
			
			

		)
	)
)

#### Server
server <- function(input, output, session) {
     rv <- reactiveValues(current_genes=c())

	# dynamic selection list: - gene to remove
	output$to_rm = renderUI({
		selectInput('to_rm2', "Select gene to remove:", c(" ", rv$current_genes), " ")
	})

	# Input -> Dotplot
	plotDotPlot <- function() {
		features <- rv$current_genes
		if (length(features) == 0) {
			features <- default_genes
		}
		pt = My_DotPlot_shiny(input_data$mean, input_data$detect, features=features, type_class=input$subset)
		print(pt)
	}


	# Make DotPlot
	output$DotPlot <- renderPlot({
		plotDotPlot()}, height=800, width=500
	)

	# Input -> UMAP
	plotUMAP <- function() {
		pt = Dim_plot(input_data, feature=input$umapfeature, type_subset=cell_types[[input$subset]], label=input$labelumap)
		print(pt)
	}
	# Make UMAP
	output$UMAPPlot <- renderPlot({
		plotUMAP()}, height=400, width=500
	)

	# Write current gene set
	output$Geneset <- renderText({paste(current_genes, collapse=", ")})

	## Buttons
	# Add
	observeEvent(input$addgene, { 
		if (input$gene != default_genes_symbol) {
			rv$current_genes <- c(rv$current_genes, input$gene)
		}
	})
	# Clear
	observeEvent(input$cleargenes, {rv$current_genes <- c()})

	# Remove
	observeEvent(input$removegene, {
		if (input$to_rm2 != " ") {
			rv$current_genes <- rv$current_genes[!rv$current_genes %in% input$to_rm2]
		}
	})


	# download DotPlot plot
	output$downloadPlot <- downloadHandler(
		filename = paste("DotPlot.png", sep="_"),
		content = function(file) {
		png(file, width=8, height=8, units="in", res=input$res)
		plotDotPlot()
		dev.off()
    })  
	
	# Stop when the browser tab is closed.
	session$onSessionEnded(stopApp)
}

#### shinyApp Call
shinyApp(ui = ui, server = server)
