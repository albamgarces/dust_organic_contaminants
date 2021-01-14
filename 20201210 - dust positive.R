## Script automatically generated on Thu Dec 10 10:43:08 2020

library(patRoon)
#####
#run to see if patRoon can find all its dependecies
#patRoon::verifyDependencies()
#install and set those not found
#options(patRoon.path.MetFragCL = "C:/Program Files/MetFrag2.4.5-CL.jar")
#options(patRoon.path.pngquant = "C:/Program Files/pngquant")
#options(patRoon.path.MetFragCompTox = "C:/Users/ASUS/Documents/CompTox_17March2019_SelectMetaData.csv")
#options(patRoon.path.MetFragPubChemLite = "C:/Users/ASUS/Documents/PubChemLite_31Oct2020_exposomics.csv")
##########
# Checking ProteoWizard... found! (directory 'C:/PROGRA~1/OPENMS~1.0/share/OpenMS/THIRDP~1/pwiz-bin')
# Checking OpenMS... found!
# Checking pngquant... found! (in 'C:/Program Files/pngquant')
# Checking SIRIUS... found!
# Checking MetFrag CL... found! (in 'C:/Program Files')
# Checking MetFrag CompTox Database... found! (in 'C:/Users/ASUS/Documents')
# Checking MetFrag PubChemLite Database... found! (in 'C:/Users/ASUS/Documents')
# Checking OpenBabel... found!
###########

library(readxl)
library(rJava)
library(plyr)

# -------------------------
# initialization
# -------------------------


#workPath <- "~/12 - Normal Trials Non Target/Indoor Dust non target/2020.12.10 - dust positive"
#setwd(workPath)

anaInfo <- generateAnalysisInfo(paths = c("mzML_pos"),
                                groups = c("Blank", "Blank", "House", "Public"),
                                blanks = c("Blank", "Blank", "Blank", "Blank"))


# -------------------------
# features
# -------------------------

## --> options(patRoon.path.OpenMS = "C:/Program Files/OpenMS-2.6.0/bin")
# directory with the OpenMS binaries

# Find all features.
# NOTE: see manual for many more options
fList <- findFeatures(anaInfo, "openms")

# Group and align features between analysis
fGroups <- groupFeatures(fList, "openms")

# Basic rule based filtering
fGroups <- filter(fGroups, preAbsMinIntensity = 10000, absMinIntensity = 1000000,
                  relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                  blankThreshold = 5, removeBlanks = TRUE,
                  retentionRange = c(50, Inf), mzRange = NULL)

### se cambia preabs y absmin para obtener menos ruido

# Filter feature groups by suspects
suspFile <- read_excel("MassBankEU_Cmpds_11042017_wMS_DTXSIDs_03102017.xlsx")

#scr <- screenSuspects(fGroups, suspFile, rtWindow = 12, mzWindow = 0.005,
                     # adduct = NULL)

#############################################################
### Error: Assertion failed. One of the following must apply:
# * checkmate::checkFactor(suspects["adduct"]):
  # * Must be of type 'factor', not 'NULL'
# * checkmate::checkCharacter(suspects["adduct"]):
  # * Must be of type 'character', not 'NULL'
#############################################################
adduct_pos <- as.vector(suspFile$`Adduct Positive Ionization`)
# unique(adduct_pos)

#scr_pos <- screenSuspects(fGroups, suspFile, rtWindow = 12, mzWindow = 0.005,
                      #adduct = "[M+H]+")

#############################################################
### Error in screenSuspects(fGroups, suspFile, rtWindow = 12, mzWindow = 0.005,  : 
        #   1 assertions failed:
        #   * Variable 'suspects': Must include the elements
        # * {name}.
#############################################################
#cambiamos el nombre de la variable y al buscar en help(screenSuspects), 
#vemos otros nombres de columna obligatorios
new_names <- c("NAME" = "name", "NEUTRAL_EXACT_MASS"="neutralMass", "FORMULA"="formula", "StdInChI"="InChi")
###mirar bien como hacerlo automático #####

names(suspFile)[2]="name"
names(suspFile)[4]="neutralMass"
names(suspFile)[5]="formula"
names(suspFile)[6]="InChI"

#############################################################
###Error in getCommandWithOptPath("obabel", "obabel") : 
#Cannot find 'obabel.exe'. Either add the correct file location to
#the PATH environment variable or set 'patRoon.path.obabel' with options().
#############################################################
#instalamos OpenBabel, rJava

scr_pos <- screenSuspects(fGroups, suspFile, rtWindow = 12, mzWindow = 0.005,
                          adduct = "[M+H]+")

###mirar como automatizarlo según adduct
### ¿reconoce los distintos adducts contenidos en la misma celda?
##+ ver como filtrar para que no repita compuestos diferentes que realmente son el mismo y aparece varias veces (columnaA, filaX = columnaB, fila Y)
###
#**** en group, hace que los groups sean únicos añadiendo el .número
#***** No funcionaba porque hay que fijarse en los nombres #
#de las columnas y que sean exactamente como los piden ***** #
###


fGroups <- groupFeaturesScreening(fGroups, scr_pos)


# -------------------------
# annotation
# -------------------------


# Retrieve MS peak lists
avgPListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)

#probar con 1, 2 el clustermz
mslists <- generateMSPeakLists(fGroups, "mzr", maxMSRtWindow = 5, precursorMzWindow = NULL,
                               avgFeatParams = avgPListParams, avgFGroupParams = avgPListParams)
# uncomment and configure for extra filtering of MS peak lists
# mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL, relMSIntThr = NULL,
#                  relMSMSIntThr = NULL, topMSPeaks = NULL, topMSMSPeaks = NULL,
#                  deIsotopeMS = FALSE, deIsotopeMSMS = FALSE)

# Calculate formula candidates
formulas <- generateFormulas(fGroups, "genform", mslists, relMzDev = 5,
                             adduct = "[M+H]+", elements = "CHNOP",
                             calculateFeatures = TRUE, featThreshold = 0.75)
## buscar todos los posibles elementos en la columna fórmulas (elements="todos")
##no prioritario

### ver cómo añadir todos los elementos que encuentre

# Find compound structure candidates
compounds <- generateCompounds(fGroups, mslists, "metfrag", method = "CL", dbRelMzDev = 5,
                               fragRelMzDev = 5, fragAbsMzDev = 0.002,
                               adduct = "[M+H]+", database = "pubchem", maxCandidatesToStop = 2500)
## maxcandidates=100 para comprobar qué pasa

compounds <- addFormulaScoring(compounds, formulas, TRUE) 

# Perform automatic generation of components
components <- generateComponents(fGroups, "ramclustr", ionization = "positive")


# -------------------------
# reporting
# -------------------------


reportCSV(fGroups, path = "report", reportFeatures = FALSE, formulas = formulas,
          compounds = compounds, compoundsNormalizeScores = "max",
          components = components)

reportPDF(fGroups, path = "report", reportFGroups = TRUE, formulas = formulas, reportFormulaSpectra = TRUE,
          compounds = compounds, compoundsNormalizeScores = "max",
          components = components, MSPeakLists = mslists)

reportHTML(fGroups, path = "report", reportPlots = c("chord", "venn", "upset", "eics", "formulas"),
           formulas = formulas, compounds = compounds, compoundsNormalizeScores = "max",
           components = components, MSPeakLists = mslists,
           selfContained = FALSE, openReport = TRUE)

