#' ---
#' title: "Bioinformatics Class 5"
#' author: "Peng Fei Huang"
#' date: "April 18, 2018"
#' output:
#'  html_document:
#'    code_folding: hide
#' ---

# Bioinformatics Class 5
# Plots 

x <- rnorm(1000,0)
summary(x)

# les see this data as a graph
boxplot(x)

# Histogram 
hist(x)



# Section 1A from lab sheet
baby <- read.table("bggn213_05_rstats/weight_chart.txt", header = TRUE)

plot(baby, type="b", pch=15, cex=1.5, lwd=2, ylim=c(0,12), lty=4, main="Baby Weight")

#Section 1B 
feat <- read.table("bggn213_05_rstats/feature_counts.txt", sep="\t", header = TRUE)
par(mar=c(5, 11, 4, 2))
barplot(feat$Count, names.arg = feat$Feature, horiz = TRUE, las=2)

# Section 2A
mfcount <- read.delim("bggn213_05_rstats/male_female_counts.txt")
barplot(mfcount$Count, names.arg = mfcount$Sample, col=rainbow(10))
mycols <- rainbow(nrow(mfcount))
barplot(mfcount$Count, col=mycols)
        
# Section 2B
updown <- read.delim("bggn213_05_rstats/up_down_expression.txt")
plot(updown$Condition1, updown$Condition2, col=updown$State)
palette(c("red","green","blue"))

# How many genes went up down and around?
table(updown$State)
unique(updown$State)

#Section 2C
read.delim("bggn213_05_rstats/color_to_value_map.r")
