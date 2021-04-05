source("include/QMC.R")


VDCSequence(10, 2)


halton.sequence <- HaltonSequence(500, 2, 3)
plot(halton.sequence$base1, halton.sequence$base2)
