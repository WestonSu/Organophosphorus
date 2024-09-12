# Script generated on Tue Dec 12 21:34:26 2023
# Edited by Wenyuan Su
library(patRoon)
# -------------------------
# initialization
# -------------------------
workPath <- "E:/patRoon/OPC"
setwd(workPath)
options(patRoon.path.pwiz = "C:/ProteoWizard") # location of ProteoWizard installation folder
options(patRoon.path.SIRIUS = "D:/Program Files/sirius") # directory with the SIRIUS binaries
patRoon::verifyDependencies()

anaInfo <- generateAnalysisInfo(paths = "E:/patRoon/OPC/analyses/raw")
# Set to FALSE to skip data pre-treatment
doDataPretreatment <- TRUE
if (doDataPretreatment)
{convertMSFiles(anaInfo = anaInfo, outPath = "analyses/mzml",from = "thermo", to = "mzML", algorithm = "pwiz", centroid = "vendor")}
anaInfo <- generateAnalysisInfo(paths = "E:/patRoon/OPC/analyses/mzml",
                                groups = c("Pooled_sample", "BK-STD"),
                                blanks = c("BK-STD", "BK-STD"),
                                norm_concs = c(NA, NA))
# -------------------------
# features
# -------------------------
# Find all features
fList <- findFeatures(anaInfo, "openms", noiseThrInt = 10000, chromSNR = 3,
                      chromFWHM = 5, mzPPM = 5, minFWHM = 1, maxFWHM = 30)
# Group and align features between analyses
fGroups <- groupFeatures(fList, "openms", rtalign = TRUE)
# Basic rule based filtering
fGroupsfil <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000, relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75, blankThreshold = 3, removeBlanks = TRUE, retentionRange = NULL, mzRange = c(100,800))
# -------------------------
# componentization
# -------------------------
# Perform automatic generation of components
components <- generateComponents(fGroupsfil, "openms", ionization = "positive")
fGroupsSel <- selectIons(fGroupsfil, components, prefAdduct = "[M+H]+", onlyMonoIso = TRUE) #[M]+/[M+H]+/[M-H]-
# -------------------------
# suspect screening
# -------------------------
#This code block can be skipped when performing nontarget analysis.
suspList <- read.csv("C:/Users/SWY/Desktop/suspectlist.csv", stringsAsFactors = FALSE)
fGroupsSusp <- screenSuspects(fGroupsSel, suspList, onlyHits =TRUE)
# -------------------------
# annotation
# -------------------------
# Retrieve MS peak lists
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
# When using data independent acquisition (DIA), set precursorMzWindow to NULL.
mslists <- generateMSPeakLists(fGroupsSusp, "mzr", maxMSRtWindow = 12, precursorMzWindow = NULL, avgFeatParams = avgMSListParams, avgFGroupParams = avgMSListParams)
# Rule based filtering of MS peak lists.
mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL, withMSMS = TRUE, relMSMSIntThr = 0.02, topMSPeaks = NULL, topMSMSPeaks = 30)
# Function to compare mzs with some tolerance for characteristic fragment ions(CFIs)
equalMZ <- function(mz1, mz2) abs(mz1 - mz2) < 0.005
mslistsF <- delete(mslists, j = function(pl, grp, ana, type) {
  if (type != "MSMS")
    return(FALSE) # Only consider MS/MS peak lists
  # the CFIs of quaternary phosphonium salts (QPSs) are used as an example
  characteristic_fragments = c(108.0122, 183.0357, 185.0513, 261.0837, 262.0909)
  num_characteristic_fragments = sum(sapply(characteristic_fragments, equalMZ, pl$mz))
  return((num_characteristic_fragments) >= 2)})
# # Calculate formula candidates
formulas <- generateFormulas(fGroupsSusp, mslists, "genform", relMzDev = 10, MSMode = "msms", batchSize = 64, maxCandidates = 1000, elements = "C[100]H[200]O[10]N[10]P[3]S[3]F[10]Cl[10]Br[10]", oc = TRUE, calculateFeatures = FALSE, timeout = 60)
compounds <- generateCompounds(fGroupsSusp, mslists, "metfrag", method = "CL", dbRelMzDev = 10, fragRelMzDev = 10, fragAbsMzDev = 0.002,database = "pubchemlite", maxCandidatesToStop = 100)
# compounds <- generateCompounds(fGroupsSusp, mslists, "sirius", relMzDev = 5, fingerIDDatabase = "pubchem", elements = "CHNOP", profile = "orbitrap", token = tokenFID)
fGroupsa <- annotateSuspects(fGroupsSusp, formulas = formulas, compounds = compounds, MSPeakLists = mslists)
compounds <- addFormulaScoring(compounds, formulas, updateScore = TRUE)
# -------------------------
# reporting
# -------------------------
# Advanced report settings can be edited in the report.yml file.
report(fGroupsa, MSPeakLists = mslists, formulas = formulas, compounds = compounds, components = NULL, TPs = NULL,  parallel = TRUE, openReport = FALSE)
