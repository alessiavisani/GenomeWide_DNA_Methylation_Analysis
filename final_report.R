
setwd("/Users/aleonmac/Desktop/UNI/DNA_RNA/final_report/Input_Data")
# loading the library and getting informations about the package minfi
library(minfi)
library(minfiData)
vignette("minfi") 

## exercise 1
#Import raw data and create an object called RGset storing the RGChannelSet object 
    
    # List the files included in the folder
    --# all the idat files for each sample either green or red colur and the csv sheet summarizing the values 
    list.files("/Users/aleonmac/Desktop/UNI/DNA_RNA/final_report/Input_Data")
    
    
    SampleSheet <- read.table("SampleSheet_Report_II.csv",sep=",",header=T)
    SampleSheet
    
    
    # Set the directory in which the raw data are stored and load the samplesheet using the function read.metharray.sheet
    baseDir <- ("/Users/aleonmac/Desktop/UNI/DNA_RNA/final_report/Input_Data")
    targets <- read.metharray.sheet(baseDir)
    #target checks that the input directory there is a sample sheet and it find it; it also modifies the sample sheet included in the input folder by adding an extra column 
    targets
    
    ?read.metharray.exp # to see if my sample sheet is constructed correctly
    RGset <- read.metharray.exp(targets = targets)
    save(RGset,file="RGset.RData") #storing the RGChannelSet
    
    RGset #investigating if the RGChannelSet was created correctly

## exercise 2 
#Create the dataframes Red and Green to store the red and green fluorescences respectively 
    Red <- data.frame(getRed(RGset))
    head(Red)
    
    Green <- data.frame(getGreen(RGset))
    head(Green)

## exercise 3
#Red/Green fluorescences for my address
    
    my_address <- "47683488"
    
    # creating Type I and Type II dataframes 
    df_I  <- data.frame(getProbeInfo(RGset, type = "I"))
    df_II <- data.frame(getProbeInfo(RGset, type = "II"))
    
    # creating the objects for my table
      if (my_address %in% df_I$AddressA | my_address %in% df_I$AddressB) {
        
        # Extract the specific row from the Type I manifest
        matched_row <- df_I[df_I$AddressA == my_address | df_I$AddressB == my_address, ]
        p_type  <- "Type I"
        p_color <- matched_row$Color 
        
      } else if (my_address %in% df_II$Address) {
        
        p_type  <- "Type II"
        p_color <- "-"
        
      }
    #creating the report table and populating it
    report_table <- data.frame(
      Sample = colnames(Red),
      `Red fluor` = as.numeric(Red[my_address, ]),
      `Green fluor` = as.numeric(Green[my_address, ]),
      Type = p_type,
      Color = p_color,
      check.names = FALSE
    )
    
    report_table
    
    # Checking  for corectness if adress exists in the Type I dataframe (AddressA or AddressB columns)
    my_address %in% df_I$AddressA | my_address %in% df_I$AddressB
    my_address %in% df_II$Address. 
    # since i only have type 2 probes for my address i'm in a comparative array two channel color also for this reason i have no record of color in the manifest file designated column
    

## exercise 4
# Create the object MSet.raw 
#we use this to extract methylated or unmethylated signals
    MSet.raw <- preprocessRaw(RGset)
    MSet.raw
    
    save(MSet.raw,file="MSet_raw.RData")
    
    #creating the objects for methylation and unmenthylation
    Meth <- as.matrix(getMeth(MSet.raw))
    Unmeth <- as.matrix(getUnmeth(MSet.raw))
    head(Meth)
    head(Unmeth)
    
    
    #When you move from an RGChannelSet (RGset) to a MethylSet (via a function like preprocessRaw()),
    #the dataset undergoes a fundamental shift in perspective: from physical chip coordinates to biological CpG coordinates.
    
    #finding the cpg_name for my address
    cpg_name <- df_II$Name[df_II$Address == my_address]
    print(cpg_name)
    
    #Grabing the values for my address in the MSet
    meth_values   <- as.numeric(Meth[cpg_name, ])
    unmeth_values <- as.numeric(Unmeth[cpg_name, ])
    
    #Grabbing the values for my address in the RGset
    red_values    <- as.numeric(Red[my_address, ])
    green_values  <- as.numeric(Green[my_address, ])
    
    # Double Check: For Type II, Green MUST equal Meth, and Red MUST equal Unmeth
    identical(meth_values, green_values)
    identical(unmeth_values, red_values)
    
    #i use function identical() and not "==" because i want a boolean output not a comparison between arrays col by col
    

## exercise 5 
#Perform the following quality checks and provide a brief comment to each step: QCplot // check the intensity of negative controls using minfi //calculate detection pValues; for each sample, how many probes have a detection p-value higher than the threshold assigned to each group?

  #  5.1 QCplot
    qc <- getQC(MSet.raw) #calculating the median 
    qc
    plotQC(qc)
  
  # 5.2 check the intensity of negative controls using minfi
    getProbeInfo(RGset, type = "Control")
    df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control"))
    
    #Address IDs strictly assigned to "NEGATIVE" controls
    neg_addresses <- df_TypeControl$Address[df_TypeControl$Type == "NEGATIVE"]
    
    #Ensure addresses exist in our raw intensity data frames
    valid_neg_addresses <- intersect(neg_addresses, rownames(Green))
    neg_green_raw       <- Green[valid_neg_addresses, ]
    neg_red_raw         <- Red[valid_neg_addresses, ]
    
    # Convert raw intensities to log2 scale to match the standard QC metric evaluation
    neg_green_log2 <- log2(neg_green_raw)
    neg_red_log2   <- log2(neg_red_raw)
    
    # Calculate the mean negative control intensity per sample for both channels
    mean_neg_green <- colMeans(neg_green_log2, na.rm = TRUE)
    mean_neg_red   <- colMeans(neg_red_log2, na.rm = TRUE)
    
    # summary table
    neg_control_summary <- data.frame(
      Mean_Log2_Green = round(mean_neg_green, 3),
      Mean_Log2_Red   = round(mean_neg_red, 3))
    
    neg_control_summary
  
## exercise 6
#calculate detection pValues; for each sample, how many probes have a detection p-value higher than the threshold assigned to each group?
    
