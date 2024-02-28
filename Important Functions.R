#Get Current Working Directory
getwd()
#Set New Working Directory
setwd("I:/Bioinformatics/Scripts")
-----------------------------------------------------------

#If you want to check what things you have done througout the full session at a particular time, you should use 
sessionInfo()
-----------------------------------------------------------

#TO SEE THE STRUCTURE OF DATA
str(datafilename)
-----------------------------------------------------------

#TO SEE THE NUMBER OF ROWS IN DATA
nrow(datafilename)
-----------------------------------------------------------

# CHECK  VERSIONS (Bioconductor, Packages, R)
R.version
version()
packageVersion("BiocManager")
packageVersion("dplyr")
# CHECK BIOCONDUCTOR VERSION
tools:::.BioC_version_associated_with_R_version()

-----------------------------------------------------------

# HOW TO GET HELP FOR ANY PACKAGES? HOW TO SEE R DOCUMENTATION FOR A PACKAGE? 
# HOW TO BROWS VIGNETTES
?pacakgename
help(packagename)
browseVignettes("packagename")
?BSgenome
help(BSgenome)
browseVignettes("BSgenome")
help(package="GenomicRanges", help_type = "html")
-----------------------------------------------------------

# WHEN YOU WANT TO KNOW ABOUT A PARTICULAR FUNCTION, JUST DO THIS FILLOWING WAY-
packagesname::functionname
-----------------------------------------------------------
# KEYBOARD SHORTCUT
ctrl + shift + M = Pipe operator (%>%)
alt + - = Assignment operator (<-)
ctrl + shift + R = Sectioncreate hobe
-----------------------------------------------------------
# Install packages from source file and by downloading the packages file
# download the source file
download.file(
"https://github.com/al2na/methylKit/releases/download/v0.99.2/methylKit_0.99.2.tar.gz",
               destfile="methylKit_0.99.2.tar.gz")
# install the package from the source file
install.packages("methylKit_0.99.2.tar.gz",
                 repos=NULL,type="source")
# delete the source file
unlink("methylKit_0.99.2.tar.gz")
-----------------------------------------------------------
library(MASS)
ls("package:MASS") # functions in the package
ls() # objects in your R enviroment
# get help on hist() function
?hist
help("hist")
# search the word "hist" in help pages
help.search("hist")
??hist
-----------------------------------------------------------
# How to save a R object for future use
# save your R objects
save(cpgi.df,enh.df,file="mydata.RData")
load("mydata.RData")

# saveRDS() can save one object at a type
saveRDS(cpgi.df,file="cpgi.rds")
x=readRDS("cpgi.rds")#When you load this rds file you have to assign it into new variable
head(x)

save(x, file = "xdata.RData")
load("xdata.RData")
x

# Dataframe
mydata <- data.frame(chr = c("chr1", "chr1", "chr2", "chr2"),
                     strand = c("-","-","+","+"),
                     start = c(200,4000,100,400),
                     end=c(250,410,200,450))
#change column names
names(mydata) <- c("chr","start","end","strand")
mydata # OR this will work too

saveRDS(mydata, file = "mydata.rds")
mydata = readRDS("mydata.rds")

save(mydata, file = "mydata.RData")
load("mydata.RData") #when we load it we need to specify the directory
mydata 
-----------------------------------------------------------
# How to delete file from directory
# At first go to the target directory where files you want to delete, then follow the below codes
unlink("filename")
unlink(c("filename", "filename", "filename"))



