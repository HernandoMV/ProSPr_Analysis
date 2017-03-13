#ProSPr_explore.R
#Explore the ProSPr resource. Please follow the instructions

#1      ####Load datasets and functions####
#You need to specify the working directory here, along with the paths of the datasets you want to analyze.
#The paths that you might need to specify are labelled with *** at the end

#working directory
setwd("/Users/vergara/Desktop/EMBL/Images/PrImR6/5PrImR6_Atlas/Feb17/ProSPr_Explore/RestOfAnimal/Analysis") #***

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

#run shiny:
shinyApp(ui = ui, server = server)

##################################PASTE SELECTION########################################################

#paste selection:
mySelection <- 
  "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349"
#generate a random string to save in a unique folder
ClusterName <- MHmakeRandomString()
dir.create(ClusterName)

#define cells to plot and supervoxels
cells2plot = rownames(Clust_Prof)[unique(as.numeric(strsplit(mySelection,',')[[1]]))]
#Get the SuperVoxels in these cells
SuperVoxel_List <- NULL
for (CellName in cells2plot){
  SuperVoxel_List = append(SuperVoxel_List,strsplit(as.character(Clust_DF[Clust_DF$clustersID==CellName,]$supervoxelsID),",")[[1]])
}


#label t-sne plot:
p <- LabelCellsInTsne(Clust_Prof,tsne_DF,cells2plot,ClusterName)
print(p)
ggsave(paste(ClusterName,"/tSNE_",ClusterName,".pdf",sep=""), width=10, height=10, dpi=700, useDingbats=FALSE)


#create text file with cell names:
fileConn<-file(paste(ClusterName,"/Cells_",ClusterName,".txt",sep=""))
writeLines(cells2plot, fileConn)
close(fileConn)


#create heatmap:
#pdf(paste(ClusterName,"/Heatmap_",ClusterName,".pdf",sep=""))
#hv <- Create_Heatmap(Clust_Prof,cells2plot,ClusterName)
#dev.off()
Create_Heatmap_Pretty(Clust_Prof,cells2plot,ClusterName,paste(ClusterName,"/Heatmap_",ClusterName,".pdf",sep=""))

#create 3d visualization:
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
DesiredGeneNames <- c("sert")
minCoexp <- 1

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
    #Get the SuperVoxels in these cells
    SuperVoxel_List <- NULL
    for (CellName in cells2plot){
      SuperVoxel_List = append(SuperVoxel_List,strsplit(as.character(Clust_DF[Clust_DF$clustersID==CellName,]$supervoxelsID),",")[[1]])
    }
    
    
    #label t-sne plot:
    p <- LabelCellsInTsne(Clust_Prof,tsne_DF,cells2plot,ClusterName)
    print(p)
    ggsave(paste(ClusterName,"/tSNE_",ClusterName,".pdf",sep=""), width=10, height=10, dpi=700, useDingbats=FALSE)
    
    
    #create text file with cell names:
    fileConn<-file(paste(ClusterName,"/Cells_",ClusterName,".txt",sep=""))
    writeLines(cells2plot, fileConn)
    close(fileConn)
    
    
    #create heatmap:
    #pdf(paste(ClusterName,"/Heatmap_",ClusterName,".pdf",sep=""))
    #hv <- Create_Heatmap(Clust_Prof,cells2plot,ClusterName)
    #dev.off()
    Filename = paste(ClusterName,"/Heatmap_",ClusterName,".pdf",sep="")
    Create_Heatmap_Pretty(Clust_Prof,cells2plot,ClusterName,Filename)
    
    
    
    #create 3d visualization:
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

#Get the SuperVoxels in these cells
SuperVoxel_List <- NULL
for (CellName in cells2plot){
  SuperVoxel_List = append(SuperVoxel_List,strsplit(as.character(Clust_DF[Clust_DF$clustersID==CellName,]$supervoxelsID),",")[[1]])
}

#label t-sne plot:
p <- LabelCellsInTsne(Clust_Prof,tsne_DF,cells2plot,ClusterName)
print(p)
ggsave(paste(dirname(filepath),"/tSNE_",ClusterName,".pdf",sep=""), width=10, height=10, dpi=700, useDingbats=FALSE)

#create heatmap:
#pdf(paste(ClusterName,"/Heatmap_",ClusterName,".pdf",sep=""))
#hv <- Create_Heatmap(Clust_Prof,cells2plot,ClusterName)
#dev.off()
Create_Heatmap_Pretty(Clust_Prof,cells2plot,ClusterName,paste(dirname(filepath),"/Heatmap_",ClusterName,".pdf",sep=""))

#create 3d visualization:
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







