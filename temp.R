library(ggplot2)

# create fictitious data
a <- runif(10)
c <- runif(7)

# data groups
group <- factor(rep(1:2, c(10, 7)))

# dataframe
mydata <- data.frame(c(a,c), group)
names(mydata) <- c("value", "group")

# function for computing mean, DS, max and min values
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# ggplot code
p1 <- ggplot(aes(y = value, x = factor(group)), data = mydata)
p1 <- p1 + stat_summary(fun.data = min.mean.sd.max, geom = "boxplot") + 
  geom_jitter(position=position_jitter(width=.2), size=3) + 
  ggtitle("Boxplot con media, 95%CI, valore min. e max.") + 
  xlab("Gruppi") + 
  ylab("Valori")


