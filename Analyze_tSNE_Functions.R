#Analyze_tSNE_Functions.R
#This script contains functions for the tSNE analysis of ProSPr

########create shiny object#########
library(shiny)
library(DT)
server <- function(input, output) {
  # render the plot
  output$plot1 <- renderPlot({
    plot(tsne_DF$V1, tsne_DF$V2)
  })
  
  # set the options for the brush technique
  output$plotui <- renderUI({
    plotOutput("plot1", height=1000,
               brush = brushOpts(id = "plot_brush")
    )
  })
  
  output$plot_brushed_points <- renderPrint({
    df <- tsne_DF
    # this function gets the data for you
    res <- brushedPoints(df, input$plot_brush, "V1","V2") 
    # mpg = name of x variable, disp = name of y variable
    paste(rownames(res), sep="", collapse=",")
    # puts the results in a datatable format
  })
}

ui <- fluidPage(
  # render the plot
  uiOutput("plotui"),  
  #dataTableOutput("plot_brushed_points")
  verbatimTextOutput("plot_brushed_points")
)

########Random String#########
#Define a function to create a random string
MHmakeRandomString <- function(n=1, lenght=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}


########Label Cells in tSNE#########
#CellMat is the matrix of expression values for each cell
#CellsToLabel is a list of the names of the cells to label
#tSNEdf is a data frame containing the position of the dots (cells) in the tSNE
#Ptit is the title of the plot
LabelCellsInTsne <- function(CellMat, tSNEdf, CellsToLabel, Ptit){
  colsforSC <- rep("gray",nrow(CellMat))
  names(colsforSC) <- rownames(CellMat)
  colsforSC[CellsToLabel] = 'red'
  p <- ggplot(tSNEdf, aes(x = V1, y = V2, color = colsforSC)) +
    scale_colour_manual(values = c("grey","red")) +
    geom_point(size = 3, alpha = 0.9) + 
    ggtitle(ClusterName) +
    theme(axis.title=element_text(size=14)) +
    theme(plot.title = element_text(size=18,vjust=1.1, face="bold")) +
    theme(axis.title=element_text(size=14)) +
    theme(legend.justification=c(0,0), legend.position='none',#legend.position=c(.9,0.73),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 14, colour = "black", angle = 0),
          legend.title = element_text(size = 0),
          legend.box = "horizontal") +
    theme(legend.key = element_blank()) +
    theme(legend.background = element_rect(fill=alpha('grey', 0.0)))
  return(p)
}



########Create heatmap of gene expression#########
library(ggplot2)
library(vegan)
library(gplots)
# set custom distance and clustering functions
hclustfunc <- function(x) hclust(x, method="complete")
#distfunc <- function(x) dist(x,method="euclidean")
#distfunc <- function(x) as.dist(1-cor(t(x),method="spearman"))
distfunc <-function(x) vegdist(x, method="jaccard")
#distfunc <-function(x) bcdist(x,rmzero = TRUE)
my_palette <- colorRampPalette(c("black","black", "yellow", "red"))(n = 95)

#CellMat is the matrix of expression values for each cell
#CellsToLabel is a list of the names of the cells to label
#Ptit is the title of the plot
Create_Heatmap <- function(CellMat, CellsToLabel, Ptit){
  submat <- CellMat[CellsToLabel,]
  submat = submat[,colSums(submat)>0]  
  hv<-heatmap.2(submat, cexCol=0.5, cexRow=0.5, trace='none', main = Ptit,
                hclust = hclustfunc, 
                distfun= distfunc,
                col=my_palette,
                margins=c(6, 6),
                #lwid controls the 3 relative sizes of the dendrogram, colors, and matrix
                #lhei the same in the vertical direction
                #lmat is something I don't quite get :)
                #RowSideColors=colsforrows,
                #lmat=rbind(c(5,0,4),c(3,1,2)),
                #lwid=c(1,5,1), 
                #lhei=c(0.15,0.9),
                dendrogram="row",
                key.title=NA
                #labRow=newrownames
  )
  return(hv)
}

########Create heatmap of gene expression with a better visualization#########
#create heatmap:
library(pheatmap)
library(vegan)

Create_Heatmap_Pretty <- function(CellMat, CellsToLabel, Ptit, NameToSave){
  submat <- CellMat[CellsToLabel,]
  submat = submat[,colSums(submat)>0]  
  hv<-pheatmap(submat, cellwidth = 2, cellheight=2,
               cluster_rows= hclust(vegdist(submat, method="jaccard"), method = "complete"),
               cluster_cols= hclust(vegdist(t(submat), method="jaccard"), method = "complete"),
               main = Ptit, fontsize =2.3, legend = T,
               filename = NameToSave
  #clustering_distance_rows="canberra", 
  #clustering_distance_cols="canberra")
  )
  dev.off()
  #return(hv)
}


########Plot SuperVoxels in 3D#########
library(rgl)
#Reference3D are the 'background' voxels
#Coordinates3D are the 3D coordinates of the supervoxels
#SVtoHighligth is a list of the SV that will be highlighted
PlotSV3D <- function(Reference3D, Coordinates3D, SVtoHighligth){
  plot3d(Reference3D,size=6,col='grey92',alpha=0.05)#,box=NA)
  spheres3d(Coordinates3D[rownames(Coordinates3D) %in% SVtoHighligth,],radius=2.18,col="red",alpha=0.6)
  aspect3d("iso")
  par3d("windowRect"= c(0,0,800,800))
  bg3d("black")
  rgl.pop(type = "bboxdeco") #remove the box
}

########Get a list of cells expressed in every gene#########
#GeneNames is a list of names that should match the column names of Exp_matrix, which contains a matrix
#with values of expression for cells (cells x genes)
GenesToCells <- function(GeneNames, Exp_matrix){
  #check that all genes are in the matrix
  if (!(all(DesiredGeneNames %in% colnames(Exp_matrix)))){
    stop('Some gene names do not match the column names of the matrix')
  }
  #initialize the list
  ListsOfCells <- vector("list", length(GeneNames))
  ListsOfCells <- setNames(ListsOfCells, GeneNames)
  #fill the list
  for (Gene in GeneNames){
    ListsOfCells[[Gene]] <- names(list(which(Exp_matrix[,Gene]>0))[[1]])
  }  
  return(ListsOfCells)
}


