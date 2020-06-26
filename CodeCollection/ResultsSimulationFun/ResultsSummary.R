
# 14.04.2020

library(ggplot2)
library(reshape2)

Rand = read.table("RandIndexResults.txt", header = T)

# Reorder variables:

Rand = Rand[ , c(1, 2, 7, 5, 3, 4, 8, 6)]

round(colMeans(Rand), 3)

round(apply(Rand, 2, sd), 3)

Data = melt(Rand)

Norm = c(rep("L22", 200), rep("L2", 200))

Data = cbind(Data, Norm)

colnames(Data) = c("Method", "ARI", "Norm")

Labs = c(expression(L[2]^2 + LU + lambda[out]), expression(L[2]^2 + LU),
         expression(L[2]^2 + lambda[out]), expression(L[2]^2),
         expression(L[2] + LU + lambda[out]), expression(L[2] + LU),
         expression(L[2] + lambda[out]), expression(L[2]))

ggplot(Data, aes(x = Method, y = ARI, fill = Norm)) + 
  geom_boxplot() +
  scale_x_discrete(labels=Labs) + 
  scale_fill_manual(values = c("#FFFFFF", "#00FFFF")) +
  theme(legend.position = "none", axis.text=element_text(size=7),
        axis.title=element_text(size=10))

ggsave('ARI.png', device = 'png', dpi = 600, width = 4.5, height = 2.5)