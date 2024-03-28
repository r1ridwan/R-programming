# Get Current Working Directory
getwd()


# Set New Working Directory
# Change working directory	(Ctrl+Shift+H)
# From "Session", click on 'set working directory' then 'choose directory'
setwd("I:/Bioinformatics/Scripts")


# To check full session at a particular time
sessionInfo()


# TO SEE THE STRUCTURE OF DATA
str(datafilename)


# TO SEE THE NUMBER OF ROWS IN DATA
nrow(datafilename)


# CHECK  VERSIONS (Bioconductor, Packages, R)
R.version
version()
packageVersion("BiocManager")
packageVersion("dplyr")
BiocManager::version()


# GET HELP/Documentation FOR ANY PACKAGE
?pacakgename
?dplyr
help(packagename)
help(dplyr)
help(package="dplyr", help_type = "html")


# HOW TO BROWS VIGNETTES
browseVignettes("packagename")
browseVignettes("BSgenome")


# TO KNOW ABOUT A PARTICULAR FUNCTION, JUST DO THIS FILLOWING WAY-
packagesname::functionname
dplyr::filter


# KEYBOARD SHORTCUTS FOR R STUDIO
# https://support.posit.co/hc/en-us/articles/200711853-Keyboard-Shortcuts-in-the-RStudio-IDE
ctrl + shift + M = Pipe operator (%>%)
alt + - = Assignment operator (<-)
ctrl + shift + R = Sectioncreate hobe
Move cursor to Console =	Ctrl+2
Clear console =	Ctrl+L	
Move cursor to beginning of line =	Home	
Move cursor to end of line =	End	
Navigate command history =	Up/Down	
Popup command history =	Ctrl+Up	
Interrupt currently executing command	= Esc	
Change working directory	= Ctrl+Shift+H


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


library(MASS)
ls("package:MASS") # functions in the package
ls() # objects in your R enviroment
# get help on hist() function
?hist
help("hist")
# search the word "hist" in help pages
help.search("hist")
??hist


## How to save a R object for future use?
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

# Create Data frame and Save it
mydata <- data.frame(chr = c("chr1", "chr1", "chr2", "chr2"),
                     strand = c("-","-","+","+"),
                     start = c(200,4000,100,400),
                     end=c(250,410,200,450))
# change column names
names(mydata) <- c("chr","start","end","strand")
mydata # OR this will work too

saveRDS(mydata, file = "mydata.rds")
mydata = readRDS("mydata.rds")

save(mydata, file = "mydata.RData")
load("mydata.RData") # when we load it we need to specify the directory
mydata 


# How to delete file from directory
# At first go to the target directory where files you want to delete, then follow the below codes
unlink("filename")
unlink(c("filename", "filename", "filename"))


