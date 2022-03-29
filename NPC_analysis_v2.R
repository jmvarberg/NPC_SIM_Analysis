#Check if the folder "graphs" exists in the current directory, if not creates it
ifelse(!dir.exists("graphs"), dir.create("graphs"), "Graphs folder exists already")
#Check if the folder "output" exists in the current directory, if not creates it
ifelse(!dir.exists("output"), dir.create("output"), "Output folder exists already")

## load required packages
if (!require(install.load)) 
{ install.packages("install.load"); library(install.load) 
}
install_load("plyr", "tidyverse", "geometry", "gtools", "data.table", "ggplot2", "reshape2", "data.table", "ggforce", "cowplot", "bio3d", "purrr", "dbscan")

#check to see if a cell cycle staging output file is in directory, and if it is use for analysis
teststage = list.files(path=".", pattern="cell_cycle_staging.csv")
if(length(teststage) > 0) {
    staging <- fread(teststage, drop="V1")
} else {
    print("No staging file found")
}

##define custom functions

#Function to convert from pixel space to micron space. Can change values to whatever correct X, Y, and Z pixel sizes are. The defaults are setup according to OMX SIM image pixel size of 40nm with Z-step size of 125n
NPC_pix2microns <-function(dataset) {
  dataset %>% 
    mutate(x = x*0.0400015,
           y = y*0.0400015,
           z = z*0.125)
}

#Function to calcualte an object's sphericity from user-defined volume and surface area values
sphericity <- function(volume, area) {
  (pi^(1/3) * (6*volume)^(2/3))/area
}

##Custom function to extract surface area from 3D convexhull with error handling. 
surf.area <- function(x, cond="problem") {
  out <- tryCatch(
    {
      convhulln(x, options="FA")$area
    }, error = function(cond) {
      message("Error with convex hull computation")
      return(NA)
    }, warning = function(cond) {
      message("Warning during convex hull computation")
      return(NA)
    }
  )
  return(out)
}

##Custom function to extract volume from 3D convexhull with error handling. 
volume <- function(x, cond="problem") {
  out <- tryCatch(
    {
      convhulln(x, options="FA")$vol
    }, error = function(cond) {
      message("Error with convex hull computation")
      return(NA)
    }, warning = function(cond) {
      message("Warning during convex hull computation")
      return(NA)
    }
  )
  return(out)
}

#Function to calculate surface area, volume, sphericity and number of NPCs for dataframe of x/y/z NPC coordinates. Returns as a data frame. Uses convhulln function from 'geometry' package.
NPC.stats <-function(x) {
  num.NPCs <- nrow(x)
  SA <- surf.area(x)
  VOL <- volume(x)
  sphericity.x <- sphericity(VOL, SA)
  density <- num.NPCs/SA
  data.frame(SA, VOL, sphericity.x, num.NPCs, density)
}

#Function to optimize 3D convex hull for surface area and volume calculations. Allows removal of 10 percent of points to help detect and remove outliers from background.
#Not sure how this handles datasets that throw convhull errors. Can try to replace convhulln calls with surf.area and volume functions but then not sure how it will handle NA values.
hullopt <- function(data) {
    input <- data
    centroid <- data %>% 
      summarise_all(mean)
    #print(centroid)
    input$dist2cent <- as.numeric(dist.xyz(as.matrix(input), as.matrix(centroid))) #dist.xyz has to have matrices as input
    input <- input %>% arrange(desc(dist2cent))
    rows <- nrow(input)
    loop <- 1
    
    while (0.1*rows > loop) {
        i <- loop+1
        current <- input[loop:rows, ]
        currenthull <- convhulln(current[,1:3])
        currentrows <- nrow(current)
        currentpct <- length(unique(c(currenthull)))/currentrows
        currentSA <- convhulln(current[,1:3], options="FA")$area
        currentV <- convhulln(current[,1:3], options="FA")$vol
        currentsph <- sphericity(currentV, currentSA)
        test <- input[i:rows, ]
        testhull <- convhulln(test[,1:3])
        testrows <- nrow(test)
        testpct <- length(unique(c(testhull)))/testrows
        testSA <- convhulln(test, options="FA")$area
        testV <- convhulln(test, options="FA")$vol
        testsph <- sphericity(currentV, currentSA)
        
        loop <- loop + 1
        #print(paste(loop-1, currentpct, testpct, currentsph, testsph))
        
        if(testpct < currentpct) {
            loop <- loop-1
            break
        }
    }
    
    output <- input[loop:rows, 1:3]
    optinfo <- loop-1
    return(list("info" = optinfo, "output" = output))
    
}

#Function to perform DBSCAN cluster analysis. input is points after hull optimization.
NPC.dbscan <- function(x) {
  d <- dbscan(x, eps=0.2, minPts=3)
  c <- as.vector(d$cluster)
  total <- length(c)
  clustd <- length(c[c>0])
  frac <- clustd/total
  data.frame(total, clustd, frac)
}


# #Function to calculate all pair-wise NPC-NPC distances, then ask for each NPC, how many NPCs are within a short distance. Output average and standard deviation for each nucleus.
# clustered <- function(data) {
#     distmat <- dist.xyz(data, data)
#     distmat[distmat==0] <- NA
#     clust <- as.data.frame(rowSums(distmat <= 0.2, na.rm=TRUE))
#     total <- nrow(clust)
#     positives <- colSums(clust > 2, na.rm=TRUE)
#     clustfract <- positives/total
#     return(clustfract)
# }


# Analysis for NPC stats - density, surface area, volume, shape

my.files <- list.files(pattern = "\\.xls$", path = "./data/", recursive = FALSE, full.names = TRUE) #make a list of all of the files to be analyzed. This will only include files ending in 
my.files.ordered <- mixedsort(my.files)
NPC.coordinates <- lapply(my.files.ordered, fread, select = c("x", "y", "z")) #read in x/y/z coordinates from batch.track macro results files
names(NPC.coordinates) <- gsub("\\.xls$", "", my.files.ordered) #remove the .xls from the names of the files
convert <- lapply(NPC.coordinates, NPC_pix2microns) #convert the coordinates from pixels into micron space using custom function 'NPC_pix2microns', in which pixel/slice size is provided
convert <- lapply(convert, hullopt) #run converted x/y/z coordinates for NPCs through hull optimization function to remove any outlier spots prior to surface area/volume calculations.
converted <-  map(convert, 2) #extract the dataframe objects containing filtered x/y/z data sets for NPC stats analysis later
removed <- map(convert, 1) #make list of dbls containing number of points removed during hull optimization to be added to final results table later
removed <- unlist(removed) #convert list to a dataframe so can cbind later

NPC.intensities <- lapply(my.files.ordered, fread, select = c("intensity")) #read in the intensity values for each NPC spot detected from the batch.track macro
names(NPC.intensities) <- gsub("\\.xls$", "", my.files.ordered) #remove the .xls from the names of the files
total.intensities <- lapply(NPC.intensities, sum) #calculate total intensities by adding all intensity values frome each NPC using function 'sum'
total.intensities <- ldply(total.intensities, data.frame) #convert the list of total intensities from a list of dataframes to a single dataframe
colnames(total.intensities) <- c(".id", "Sum Intensity") #rename the columns of resulting dataframe

results <- lapply(converted, NPC.stats) #calculate NPC number, NE surface area/volume, NPC density, sphericity using the custom 'NPC.stats' function on the data converted to microns
results.combined <- ldply(results, data.frame) #convert list of NPC stats to a single data frame with one row per file/image analyzed
npcdb <- lapply(converted, NPC.dbscan) #perform DBSCAN clustering analysis on hull opt xyz points
npcdb <- ldply(npcdb, data.frame) #convert DBSCAN output to a single dataframe object
colnames(npcdb) <- c("ID", "Total.NPCs", "Clust.NPCs", "Frac.Clust")

final.results <- cbind(results.combined, total.intensities$`Sum Intensity`)
final.results$Points.Removed <- removed #add number of points removed column
final.results$Clust.NPCs <- npcdb$Clust.NPCs
final.results$Frac.Clust <- npcdb$Frac.Clust
colnames(final.results) <- c("Image.Analyzed", "Surface.Area", "Volume", "Sphericity", "Number.NPCS", "NPC.Density", "Total.Intensity", "Points.Removed", "Clust.NPCs", "Frac.Clust") #change the names of the columns in the results file
final.results$Avg.Intensity <- final.results$Total.Intensity/final.results$Number.NPCS #add a new column calculating the average intensity per NPC spot
final.results$Image.Analyzed <- substr(final.results$Image.Analyzed, 9, nchar(final.results$Image.Analyzed))

#Calculate and add proximity ratio/clustering analysis values to each nucleus
my.pcfs <- list.files(pattern = "\\pcf.txt$|\\pcf.csv$", path = "./data/", recursive = FALSE, full.names = TRUE)
my.pcfs <- mixedsort(my.pcfs)
PCFs <- lapply(my.pcfs, fread)

if(length(PCFs) !=0) {

  #Convert x-values in PCF lists to nm's and calculate Proximity Ratio.
  proxRatio <- function(df) {
    colnames(df) <- c("x", "y")
    df <- df %>% 
      mutate(x = x*0.040)
    
    width <- max(subset(df, y > 0))/2 #calculate "radius", 1/2 of the max distance between two points, interpreted by non-zero values of r in PCF data
    mid_start <- (width*pi/6) - 0.16 ##calculate start as radius# pi/6, then subtract 160nm to start window before the point of divergence
    mid_stop <- mid_start+0.16 #based on reported divergence of ripley's K function between spherical point data analyzed using 2D or 3D approaches, sets window ~0.52 x radius
    roi2_range <- paste(mid_start, mid_stop, sep = "-")
    
    roi1.sum <- df %>% 
      filter(x >= 0.160 & x <= 0.320) %>%
      summarize(total = sum(y))
    
    roi2.sum <- df %>%
      filter(x >=mid_start & x <= mid_stop) %>%
      summarize(total = sum(y))
    
    ratio <- roi1.sum/roi2.sum
    
    data.frame(roi1.sum, roi2.sum, ratio, roi2_range)
    
  }
  
  PCFs.proxRatio <- lapply(PCFs, proxRatio)
  names(PCFs.proxRatio) <- my.pcfs
  PCFs.proxRatio <- ldply(PCFs.proxRatio, data.frame)
  colnames(PCFs.proxRatio) <- c("Image.Analyzed","160-320nm", "Middle", "Prox.Ratio", "Middle_Range")
  PCFs.proxRatio <- PCFs.proxRatio %>%
    mutate(Image.Analyzed = str_replace(Image.Analyzed, regex("(_pcf).*", ignore_case=TRUE), "")) %>%
    mutate(Image.Analyzed = str_replace(Image.Analyzed, "./data//", "")) %>%
    mutate_if(is.numeric, ~round(., 3)) %>%
    filter(!Prox.Ratio == "Inf")
  #left join to merge PCF info in for all nuclei in the comb.df.staged dataframe.
  final.results.pcf<- left_join(final.results, PCFs.proxRatio, by = "Image.Analyzed")
  
  write.csv(final.results.pcf, file = "./output/SIM_NPC_Quantitation_Analysis_output.csv") #write out results file as CSV to working directory
}

if(length(PCFs) == 0) {
  write.csv(final.results, file = "./output/SIM_NPC_Quantitation_Analysis_output.csv") #write out results file as CSV to working directory
}

#Add staging information if present
if(exists("staging")) {
  staging <- staging %>%
    rename(Image.Analyzed = Label)
  
  final.results.staged <- left_join(final.results.pcf, staging, by="Image.Analyzed")
  ##final.results.staged[final.results.staged == 0] <- NA
  
  #Add information about cell cycle stage to results dataframe. Also filters to remove nuclei that
  #weren't able to be staged due to cell going out of FOV, or due to defective segmentation.
  final.results.staged <- final.results.staged %>%
    filter(Length > 5) %>%
    mutate(Stage = case_when(
      (Length < 9.5 & Septated == 0 & Binucleate == 0) ~ "Early G2",
      (Length >= 9.5 & Length < 11 & Septated == 0 & Binucleate == 0) ~ "Mid G2",
      (Length >= 11  & Septated == 0 & Binucleate == 0) ~ "Late G2/Early M",
      (Length >= 11 & Septated == 0 & Binucleate == "yes") ~ "Late Mitosis",
      (Septated == "yes") ~ "G1/S")) %>% 
    mutate_if(is.numeric, round, 3)
  
  final.results.staged$Stage <- factor(final.results.staged$Stage, levels = c("Early G2", "Mid G2", "Late G2/Early M", "Late Mitosis", "G1/S"))
  #save out final results as csv
  write.csv(final.results.staged, "./output/Full_NPC_Analysis_Staged_output.csv")
  
  #save out final results filtered for mid G2 for comparisons between conditions/strains.
  
  midg2results <- final.results.staged %>% 
    filter(Stage == "Mid G2")
  write.csv(midg2results, "./output/MidG2_Results.csv")
  
  #Generate files and plots of summary stats by group
  
  counts <- final.results.staged %>% 
    group_by(Stage) %>% 
    count()
  stats <- final.results.staged %>% 
    group_by(Stage) %>% 
    summarise_if(is.numeric, list(mean = mean, sd = sd), na.rm=TRUE)
  
  sum.stats <- left_join(counts, stats, by = "Stage") %>% 
    mutate_if(is.numeric, round, 3)
  write.csv(sum.stats, "./output/Summary_statistics_by_stage.csv")
  
  #Generate Plots by group

  #NPC number
  ggplot(final.results.staged, aes(x=Stage, y=Number.NPCS)) +
    geom_sina() +
    theme_cowplot()
  ggsave2("./output/NPC_number_by_stage.pdf", width=8, height=8, units="in", dpi=300)
  
  #NPC density
  ggplot(final.results.staged, aes(x=Stage, y=NPC.Density)) +
    geom_sina() +
    theme_cowplot()
  ggsave2("./output/NPC_density_by_stage.pdf", width=8, height=8, units="in", dpi=300)
  
  #Proximity Ratio
  ggplot(final.results.staged, aes(x=Stage, y=Prox.Ratio)) +
    geom_sina() +
    theme_cowplot()
  ggsave2("./output/Proximity_Ratio_by_stage.pdf", width=8, height=8, units="in", dpi=300)
  
}




