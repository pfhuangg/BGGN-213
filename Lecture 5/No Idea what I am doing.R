library(bio3d)

x <- "4AKE"

s1 <- read.pdb(x)
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
plotb3(s1.b, typ="l", ylab="Bfactor") 




library(bio3d)
x <- "1AKE"
s1 <- read.pdb(x)
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
plotb3(s1.chainA$atom$b, sse=s1.chainA, typ="l", ylab="Bfactor") 


Protein_Drug <- function(x) {
  library(bio3d)
  bring_up <- read.pdb(x)
  protein.chainA <- trim.pdb(bring_up, chain="A", elety="CA")
  Graph <- plotb3(protein.chainA$atom$b, sse=s1.chainA, typ="l", ylab="Bfactor")
}


Protein <- function(x){
  library(bio3d)
  s1 <- read.pdb(x)
  s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
  plotb3(s1.chainA$atom$b, sse=s1.chainA, typ="l", ylab="Bfactor") 
}

Protein("4AKE") 
  



