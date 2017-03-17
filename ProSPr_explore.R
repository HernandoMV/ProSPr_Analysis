#ProSPr_explore.R
#Explore the ProSPr resource. Please follow the instructions

#1      ####Load datasets and functions####
#You need to specify the working directory here, along with the paths of the datasets you want to analyze.
#The paths that you might need to specify are labelled with *** at the end

#working directory
setwd("/Users/vergara/Desktop/EMBL/Images/PrImR6/5PrImR6_Atlas/Feb17/ProSPr_Explore/Just_TFs/EntireAnimal/Analysis") #***

#SPECIFIC DATASETS
#load t-sne data:
x <- load("../Cells_tSNE") #***
tsne_DF <- get(x)
rm(x)

#load clusters data frame (this contains information about the SV in the clusters or cells):
x <- load("../Cells_properties") #***
Clust_DF <- get(x)
rm(x)

#load clusters profile:
x <- load("../Cells_profile") #***
Clust_Prof <- get(x)
rm(x)


#GENERAL DATASETS
#load 3D coordinates:
load("~/Desktop/EMBL/Images/PrImR6/5PrImR6_Atlas/Nov16/SuperVoxels_3D_coordinates") #***

#get the surface points:
load("~/Desktop/EMBL/Images/PrImR6/5PrImR6_Atlas/Nov16/surfacePoints3D_6pdf") #***

#load specific functions:
source("~/Desktop/EMBL/Images/PrImR6/5PrImR6_Atlas/Scripts/ProSPr_Analysis/ProSPr_Analysis_Functions.R") #***




#2      ####Select one of the options and run it####

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ OPTION I: SELECT A SUBSET OF CELLS IN THE TSNE $$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

##################################RUN APPLICATION########################################################
#generate a random string to save in a unique folder
ClusterName <- MHmakeRandomString()
ClusterName <- MHmakeRandomString()
print(ClusterName)
#run shiny:
shinyApp(ui = ui, server = server)

##################################PASTE SELECTION########################################################
#paste selection:
mySelection <- 
  "534,645,653,654,724,911,2270,2278,2285,2286,2299,2404,2677,2805,2806,2807,2891,2894,2950,3105,3106,3180,3182,3205,3206,3207,3218,3269,725,726,913,914,926,927,937,966,985,992,1012,1082,2287,2303,2304,2418,2575,2952,3105,3106,3182,3205,3206"
#Define the number of clusters you want to have automatically in the heatmap and the 3D visualization (use 0 for default)
NumberOfClusters = 6

dir.create(ClusterName)

#define cells to plot
cells2plot = rownames(Clust_Prof)[unique(as.numeric(strsplit(mySelection,',')[[1]]))]

#label t-sne plot:
p <- LabelCellsInTsne(Clust_Prof,tsne_DF,cells2plot,ClusterName)
print(p)
ggsave(paste(ClusterName,"/tSNE_",ClusterName,".pdf",sep=""), width=10, height=10, dpi=700, useDingbats=FALSE)

#create text file with cell names:
fileConn<-file(paste(ClusterName,"/Cells_",ClusterName,".txt",sep=""))
writeLines(cells2plot, fileConn)
close(fileConn)

#create heatmap:
subClust <- Create_Heatmap_Pretty(Clust_Prof,cells2plot,ClusterName,
                      paste(ClusterName,"/Heatmap_",ClusterName,".pdf",sep=""),
                      NumberOfClusters)

#create 3d visualization:
#find the supervoxels composing the cells in each cluster
SuperVoxel_List <- getListOfSV(cells2plot,subClust,Clust_DF)
#plot them with the same colors as the heatmap
PlotSV3D(surf3d, coord3d, SuperVoxel_List)
pic_count = 0

########Save as many pictures as you want#######
pic_count = pic_count + 1
rgl.snapshot(filename = paste(ClusterName,"/3D_",pic_count,"_",ClusterName,".png",sep=""))


#########################################################################################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ OPTION II: SEARCH BY GENE(S) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#Write the name of the gene(s) here, the minimum coexp value, and run the block
DesiredGeneNames <- c("sert","Hb9")
minCoexp <- 2
#Define the number of clusters you want to have automatically in the heatmap and the 3D visualization (use 0 for default)
NumberOfClusters = 2

#control checks:
if (all(DesiredGeneNames %in% colnames(Clust_Prof))){
  #generate a random string to save in a unique folder
  ClusterName <- paste(paste(DesiredGeneNames,collapse="_"),paste('minCoexp',minCoexp,sep="",collapse=""),sep="-")
  dir.create(ClusterName,showWarnings = FALSE)
  #define cells to plot and supervoxels
  #use a function to get a list of the cells expressing every gene
  CellList <- GenesToCells(DesiredGeneNames, Clust_Prof)
  #get their intersection
  #cells2plot <- Reduce(intersect,CellList)
  #Instead of the full co-expression (above), make it more flexible to select the minimum gene coexpression
  #make a table with the counts per cell
  CellTable <- as.matrix(table(as.vector(unlist(CellList))))
  #adjust the minCoexp value if it is bigger than the list size
  if (minCoexp > length(DesiredGeneNames)){
    minCoexp = length(DesiredGeneNames)
  }
  cells2plot <- rownames(CellTable)[which(CellTable >= minCoexp)]
 
  #make a venn diagram of the coexpression overlap in the genelist is bigger than 1
  if (length(DesiredGeneNames)>1){
    library(VennDiagram)
    venn.diagram(CellList, file = paste(ClusterName,"/VennDiagram-",ClusterName,".tif",sep=""), height = 500, width = 500, margin=0.1, resolution=50, lwd = 4, cat.cex = 3, cex = 2)
  }
  #check that cells2plot is not empty, otherwise return an error, indicating that a venn diagram has been created
  if (length(cells2plot) == 0){
    stop('Sorry, no cell coexpresses these set of genes... :(   Please check the Venn-Diagram')
  }else{   
    #label t-sne plot:
    p <- LabelCellsInTsne(Clust_Prof,tsne_DF,cells2plot,ClusterName)
    print(p)
    ggsave(paste(ClusterName,"/tSNE_",ClusterName,".pdf",sep=""), width=10, height=10, dpi=700, useDingbats=FALSE)
        
    #create text file with cell names:
    fileConn<-file(paste(ClusterName,"/Cells_",ClusterName,".txt",sep=""))
    writeLines(cells2plot, fileConn)
    close(fileConn)
    
    #create heatmap:
    subClust <- Create_Heatmap_Pretty(Clust_Prof,cells2plot,ClusterName,
                                      paste(ClusterName,"/Heatmap_",ClusterName,".pdf",sep=""),
                                      NumberOfClusters)
    
    #create 3d visualization:
    #find the supervoxels composing the cells in each cluster
    SuperVoxel_List <- getListOfSV(cells2plot,subClust,Clust_DF)
    #plot them with the same colors as the heatmap
    PlotSV3D(surf3d, coord3d, SuperVoxel_List)
    pic_count = 0
  }
}else{
  #Prompt a message telling which genes are not in ProSPr
  GenesNotIn = paste(DesiredGeneNames[- which((DesiredGeneNames %in% colnames(Clust_Prof)))],collapse=',')
  stop(paste('These genes are not in ProSPr:', GenesNotIn, '; Type: "colnames(Clust_Prof)" to get all names', collapse=""))
}

########Save as many pictures as you want#######
pic_count = pic_count + 1
rgl.snapshot(filename = paste(ClusterName,"/3D_",pic_count,"_",ClusterName,".png",sep=""))



#########################################################################################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ OPTION III: ANALYZE A LIST OF CELLS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#Define the number of clusters you want to have automatically in the heatmap and the 3D visualization (use 0 for default)
NumberOfClusters = 2
#This option allows you to load a text file containing cell names, see where they are in the animal
#and in the tSNE, and cluster them separately using hierarchical clustering
#!!!!!First, create a folder and put inside that folder a text file with the cell names, each name in one
#row, then you can run this part.

#load the text file containing a set of cells
filepath <- file.choose("Choose the text file containing the cell names")
#set the cluster name according to the directory
ClusterName <- tail(unlist(strsplit(dirname(filepath), .Platform$file.sep)),1)
#open the file and get the names as a vector
filedata <- read.table(filepath)
cells2plot <- as.vector(filedata$V1)

# #Get the SuperVoxels in these cells
# SuperVoxel_List <- NULL
# for (CellName in cells2plot){
#   SuperVoxel_List = append(SuperVoxel_List,strsplit(as.character(Clust_DF[Clust_DF$clustersID==CellName,]$supervoxelsID),",")[[1]])
# }

#label t-sne plot:
p <- LabelCellsInTsne(Clust_Prof,tsne_DF,cells2plot,ClusterName)
print(p)
ggsave(paste(dirname(filepath),"/tSNE_",ClusterName,".pdf",sep=""), width=10, height=10, dpi=700, useDingbats=FALSE)

#create heatmap:
subClust <- Create_Heatmap_Pretty(Clust_Prof,cells2plot,ClusterName,
                                  paste(dirname(filepath),"/Heatmap_",ClusterName,".pdf",sep=""),
                                  NumberOfClusters)

#create 3d visualization:
#find the supervoxels composing the cells in each cluster
SuperVoxel_List <- getListOfSV(cells2plot,subClust,Clust_DF)
#plot them with the same colors as the heatmap
PlotSV3D(surf3d, coord3d, SuperVoxel_List)
pic_count = 0

########Save as many pictures as you want#######
pic_count = pic_count + 1
rgl.snapshot(filename = paste(dirname(filepath),"/3D_",pic_count,"_",ClusterName,".png",sep=""))




#########################################################################################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ OPTION IV: CREATE SUB-TNE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#This option generates a tSNE using a set of cell names, entered as a text file, and creates the datasets
#needed for its analysis using this same script with the options above

#!!!!!First, create a folder and put inside that folder a text file with the cell names, each name in one
#row, then you can run this part. Once this is run, variable names will be reassigned, so you will need
#to restart the script again

#load the text file containing a set of cells
filepath <- file.choose("Choose the text file containing the cell names")
#set the working directory
setwd(dirname(filepath))
#open the file and get the names as a vector
filedata <- read.table(filepath)
cells_subset <- as.vector(filedata$V1)
#select the subset profile
Clust_Prof_subset <- Clust_Prof[cells_subset,]
#remove genes with no information
Clust_Prof_subset = Clust_Prof_subset[,colSums(Clust_Prof_subset)>0]
#get those cells that are unique
Clust_Prof_subset <- unique(Clust_Prof_subset) #tSNE only doesn't work with duplicates
#save dataset
save(Clust_Prof_subset, file = "Clusters_subset_prof_unique")
#select the subset dataframe (contains information about the supervoxels composing the cells)
Clust_DF_subset <- Clust_DF[Clust_DF$clustersID %in% cells_subset,]
#save dataset
save(Clust_DF_subset, file = "Clusters_subset_DF")


#tSNE analysis
library(Rtsne)
xx<-vegdist((Clust_Prof_subset), method="jaccard")
set.seed(42) # Set a seed if you want reproducible results
tsne_object <- Rtsne(as.matrix(xx),initial_dims = 10, perplexity = 31, theta=0.1)  ## this is the default setting.
tsne_DF_subset <- as.data.frame(tsne_object$Y)
#save object
save(tsne_DF_subset, file="Subset_tsneDF")


#TSNE plotting the number of genes in
p <- ggplot(tsne_DF_subset, aes(x = V1, y = V2, color = apply(Clust_Prof_subset,1, function(x) length(which(x>0))))) +
  geom_point(size = 3, alpha = 0.5) + 
  scale_color_gradient(low = 'grey50', high = 'green') + 
  #scale_color_gradientn(colours = (terrain.colors(100))) + 
  #scale_colour_manual(values = c("yellow","blue","green","magenta","red")) +
  ggtitle("TSNE on Head SuperVoxels Clusters") +
  theme(axis.title=element_text(size=14)) +
  theme(plot.title = element_text(size=18,vjust=1.1, face="bold")) +
  theme(axis.title=element_text(size=14)) +
  theme(legend.justification=c(0,0), legend.position=c(.87,0.7),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 14, colour = "black", angle = 0),
        legend.title = element_text(size = 0),
        legend.box = "horizontal") +
  theme(legend.key = element_blank()) +
  theme(legend.background = element_rect(fill=alpha('grey', 0.0)))
print(p)
ggsave("TSNE_GeneRepresentation.pdf", width=10, height=10, dpi=700, useDingbats=FALSE)

#Create graphs for genes
dir.create("TSNE_GeneExpression/")
for (gene_name in colnames(Clust_Prof_subset)){
  print(gene_name)
  p <- ggplot(tsne_DF_subset, aes(x = V1, y = V2, color = Clust_Prof_subset[,gene_name]*100)) +
    geom_point(size = 3, alpha = 0.9) + 
    #scale_colour_gradientn(colours = topo.colors(5)) + 
    scale_color_gradient(low = 'grey50', high = 'red') + 
    ggtitle(gene_name) +
    theme(axis.title=element_text(size=14)) +
    theme(plot.title = element_text(size=18,vjust=1.1, face="bold")) +
    theme(axis.title=element_text(size=14)) +
    theme(legend.justification=c(0,0), legend.position=c(.87,0.7),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 14, colour = "black", angle = 0),
          legend.title = element_text(size = 0),
          legend.box = "horizontal") +
    theme(legend.key = element_blank()) +
    theme(legend.background = element_rect(fill=alpha('grey', 0.0)))
  print(p)
  ggsave(paste("TSNE_GeneExpression/TSNE_GeneExpression_",gene_name,".pdf"), width=10, height=10, dpi=700, useDingbats=FALSE)
}







