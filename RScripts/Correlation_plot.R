# The script can be used to generate a correlation matrix from data frame and plot it
# In our study, we generated a dataframe with the small intestinal signature taxa in NCGS and IBS patients and Scores of Symptoms collected from NCGS patients

# Load library
library(corrplot)

# Set working directory
setwd("path")

# load dataframe containing taxa and symptom score
my_data <- read.delim("ncgs_top_symptoms_duo.txt")

M = cor(my_data)

corrplot(M, p.mat = testRes$p, order = 'hclust', addrect = 2, sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2, pch.col = "green4", insig = 'label_sig', col = COL2('PuOr', 10))
corrplot(M, method = 'number')
corrplot(M, method = 'color', order = 'alphabet')
corrplot(M)
corrplot(M, order = 'AOE')
corrplot(M, method = 'shade', order = 'AOE', diag = FALSE)
corrplot(M, p.mat = testRes$p, method = 'square', order = 'FPC', type = 'lower', addrect = 2, sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2, pch.col = "red", insig = 'label_sig', diag = FALSE, col = COL2('PuOr', 10))
corrplot(M, method = 'ellipse', order = 'AOE', type = 'upper')
corrplot(M, order = 'hclust', addrect = 2, col = COL2('PuOr', 10))

corrplot(M, order = 'hclust', addrect = 4, col = c('#3399FF', 'orange'))

################plotting using corrplot#################
#Transform data to matrix
matrix_cor<-as.matrix(M[,-1])

corrplot(matrix_cor, outline = TRUE, tl.cex = 1.5, tl.col = "black", cl.cex = 1.5, is.corr=FALSE,
         method = "circle")

# to calculate significance
testRes = cor.mtest(Q, conf.level = 0.95)
