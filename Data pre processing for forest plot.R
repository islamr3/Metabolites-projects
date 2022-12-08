
#### Data Pre-Processing ###

library(ggplot2)
library(reshape2)
Model <- lm(cyl~., data = mtcars)
Model
coef <- Model$coefficients
coef <- melt(coef)
coef
Col1 <- rownames(coef)
Col2 <- coef[,1]
dt <- cbind(Col1, Col2)
class(dt)
dt <- as.data.frame(dt)

####
#Please try to create forest plot in r.
#
