##### Meta-Analysis ####
##사용했다는 package 들을 정리해보자
library(devtools)

###Data retrieval, OTU picking

BiocManager::install("ShortRead")

###Diversity Analysis라고 되어있는 부분

BiocManager::install("vsearch") #install.package로 하면 다른 version이라 설치가 안됨 
BiocManager::install("philr") #install.package로 하면 다른 version이라 설치가 안됨 
#                             얘 먼저 설치해야 MicrobeR 이 설치 됨
#**BiocManager::install("PhILR") #install.package로 하면 다른 version이라 설치가 안됨 **
#*package ‘PhILR’ is not available for Bioconductor version '3.14'
install_github("jbisanz/MicrobeR") #
install.packages("picante")
install.packages("zCompositions")
install.packages("lmerTest")

###Random Forest Classifiers

install.packages("randomForest")




##불러오기
library("ShortRead")
library("MicrobeR")
library("vegan")
library("picante")
library("phyloseq")
library("zCompositions")
library("base")
##library("PhILR")
library("ape")
library("lmerTest")
library("randomForest")


##### LET'S GET STARTED #####










