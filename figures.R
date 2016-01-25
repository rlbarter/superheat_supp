# code to generate the figures
library(superheat)
library(glmnet)
library(ggplot2)

source("../Manuscript/0-clean.R")




############# Generate Figure 3: Communities and crimes EDA #############

# remove extreme observations
sd <- apply(crime.df, 2, sd)
mean <- apply(crime.df, 2, mean)
upper <- mean + 3*sd
lower <- mean - 3*sd

above <- mapply(function(x, y) x > y, crime.df, upper)
above <- lapply(data.frame(above), function(x) which(x > 0))
below <- mapply(function(x, y) x < y, crime.df, lower)
below <- lapply(data.frame(below), function(x) which(x > 0))
to_rm <- unique(c(unlist(above), unlist(below)))



# scale the variables
crime.df.new <- crime.df[-to_rm,]
crime.scaled <- scale(crime.df.new)

# calculate the correaltion of each variable with the response (ViolentCrimesPerPop)
cor.vc <- cor(crime[,c("ViolentCrimesPerPop",colnames(crime.df))])[,"ViolentCrimesPerPop"]


# set the seed!
set.seed(12345681)

# perform clustering before generating heatmap
membership <- factor(kmeans(scale(crime.df), 3)$cluster)

# reverse the order of the cluster names
levels(membership) <- c(3,2,1)
# remove the extreme observations
membership.new <- as.numeric(as.vector(membership))[-to_rm]

# cluster the variables
vars.membership <- factor(kmeans(t(scale(crime.df)), 5)$cluster)
# change order of cluster names
levels(vars.membership) <- c(2,5,4,1,3)
vars.membership.new <- as.numeric(as.vector(vars.membership))

# generate the EDA heatmap for the crime dataset
crime_heat <- superheat(X = t(crime.scaled),
                       membership.cols = membership.new,
                       membership.rows = vars.membership.new,
                       heat.pal = c("white", "lightskyblue1", "navyblue", "black"),
                       heat.pal.values = c(0, 0.24, 0.3, 1),
                       order.cols = order(t(crime.scaled)["PctKidsBornNeverMar",]),
                       order.rows = order(cor.vc[-1]),
                       cluster.box = T,
                       box.size = 1,
                       box.col = "white",
                       line.size = 1,

                       # left labels
                       left.text.angle = 0,
                       left.text.size = 3,
                       left.label.size = 0.3,
                       left.heat.label = "variable",
                       left.label.text.col = rep(c("black","black"), 16),
                       left.label.pal = c("gray76","gray61"),

                       # bottom labels
                       bottom.label.text.col = c("black","white", "black"),

                       # top plot
                       yt = y[-to_rm],
                       yt.num.ticks = 4,
                       yt.plot.type = "scattersmooth",
                       yt.axis.name = "violent crimes per pop",
                       yt.point.alpha = 0.2,

                       # right plot
                       yr = cor.vc[-1],
                       yr.plot.type = "bar",
                       yr.axis.name = "Correlation with violent crime",
                       yr.num.ticks = 4,

                       # legend
                       legend.size = 4,

                       # change to print.plot = TRUE to show the plot
                       print.plot = FALSE)











######## Generate Figure 4: Communities and crimes model assessment #########

# add the response to the data frame
crime.df$ViolentCrimesPerPop <- y
# scale the data frame
crime.df <- data.frame(scale(crime.df))

# fit a lasso linear model (find regularization parameter using CV)
cv.lasso <- cv.glmnet(x = as.matrix(crime.df[,selected.vars]),
                      y = crime.df$ViolentCrimesPerPop)
lambda <- cv.lasso$lambda.min
lasso <- glmnet(x = as.matrix(crime.df[,selected.vars]),
                y = crime.df$ViolentCrimesPerPop)
# identify which variables had nonzero coefficient
selected.lasso <- which(lasso$beta[,which.min(abs(lasso$lambda - lambda))] != 0)
# refit model using only selected variables
model.lasso <- lm(ViolentCrimesPerPop ~.,
                  crime.df[,c("ViolentCrimesPerPop",names(selected.lasso))])






# log transformation - fit lasso linear model again
crime.df$ViolentCrimesPerPopLog <- log(1 + crime.df$ViolentCrimesPerPop)
cv.log.lasso <- cv.glmnet(x = as.matrix(crime.df[,selected.vars]),
                          y = crime.df$ViolentCrimesPerPopLog)
lambda.log <- cv.log.lasso$lambda.min
lasso.log <- glmnet(x = as.matrix(crime.df[,selected.vars]),
                    y = crime.df$ViolentCrimesPerPopLog)
# lasso sets only a few coefficients to zero
selected.log.lasso <- which(lasso.log$beta[,which.min(abs(lasso.log$lambda - lambda.log))] != 0)
model.log.lasso <- lm(ViolentCrimesPerPopLog ~.,
                      data = crime.df[,c("ViolentCrimesPerPopLog",names(selected.log.lasso))])




# use same membership vector as last example
membership <- as.numeric(as.vector(membership))

# select a sample of size 500
index <- sample(1:length(membership), 500)

cor.crime <- cor(t(crime.df[index,selected.lasso]))

crime_fit_heat <- superheat(X = cor.crime,
                      membership.cols = membership[index],
                      membership.rows = membership[index],
                      heat.pal = c("#F35B5B", "white", "#0D7F7F"),
                      heat.pal.values = c(0, 0.45, 1),
                      order.cols = order(model.lasso$fitted.values[index]),
                      order.rows = order(model.log.lasso$fitted.values[index]),
                      cluster.box = T,
                      box.size = 0.7,
                      line.size = 1,

                      # top plot
                      yt = model.lasso$residuals[index],
                      yt.num.ticks = 4,
                      yt.plot.type = "scatter",
                      yt.axis.name = "residuals from standard\n linear model",
                      yt.point.alpha = 0.7,

                      # right plot
                      yr = model.log.lasso$residuals[index],
                      yr.num.ticks = 4,
                      yr.axis.name = "residuals from log-transformed\n linear model",

                      # left labels
                      left.text.angle = 0,
                      left.label.text.col = c("black","white","black"),

                      # bottom labels
                      bottom.label.text.col = c("black","white","black"),

                      # legend
                      legend.size = 4,

                      # change to print.plot = TRUE to show the plot
                      print.plot = FALSE)







######################### Generate Figure 6: fMRI EDA #################################


load("fMRIdata_clean.RData")


set.seed(123)
X <- fit_feat
Y <- resp_dat


# X: images are observations, features are columns
# y: the response to each image for a single voxel

col.membership <- kmeans(Y, centers = 2)$clust
col.membership <- -(col.membership - 3)
col.membership <- factor(col.membership)
levels(col.membership) <- c("image cluster 1","image cluster 2")

row.membership <- kmeans(t(X), centers = 2)$clust
row.membership <- -(row.membership - 3)
row.membership <- factor(row.membership)
levels(row.membership) <- c("feature cluster 1","feature cluster 2")




# replace values with their ranks:
X.rank <- matrix(rank(X)/(nrow(X)*ncol(X)), ncol = ncol(X))
fmri.design <- superheat(X = t(X.rank),
                           membership.cols = col.membership,
                           membership.rows = row.membership,
                           cluster.box = T,
                           heat.pal = c("#b2182b",
                             "#ef8a62",
                             "#fddbc7",
                             "#d1e5f0",
                             "#67a9cf",
                             "#2166ac"),
                           heat.pal.values = c(0, 0.2,  0.495, 0.505, 0.8, 1),

                           # top plot
                           yt = -pc$scores[,1],
                           yt.axis.name = "First principal component\nof voxel responses",
                           yt.plot.type = "boxplot",
                           yt.plot.size = 0.7,
                           yt.point.size = 3,
                           #yt.pal = c("#ef8a62", "#67a9cf"),
                           yt.axis.size = 15,
                           yt.axis.name.size = 15,
                           yt.point.alpha = 0.6,

                           # left labels
                           left.text.angle = 0,
                           left.label.size = 0.3,
                           left.label.text.col = c("black","white"),

                           # bottom labels
                           bottom.label.text.col = c("black","white"),

                           # change to print.plot = TRUE to show the plot
                           print.plot = FALSE)



################################## Generate Figure 8: fMRI model assessment ####################################

load("fMRI_results.RData")

set.seed(123)
test.index <- sample(1:1750, 300)
row.membership <- col.membership

col.membership <- kmeans(t(Y), centers = 2)$clust
col.membership <- -(col.membership - 3)
col.membership <- factor(col.membership)
levels(col.membership) <- c("voxel cluster 1","voxel cluster 2")


colnames(Y.test) <- paste0("V", 1:ncol(Y.test))




# do the visualizations (fit new model and see correlation) on validation set
fmri.fit <- superheat(X = Y.test,
                        membership.cols = col.membership,
                        membership.rows = row.membership[test.index],
                        cluster.box = T,
                        heat.pal = c("#ffffcc",
                          "#c7e9b4",
                          "#7fcdbb",
                          "#41b6c4",
                          "#2c7fb8",
                          "#253494"),
                        heat.pal.values = c(0, 0.45,  0.495, 0.505, 0.55, 1),
                        order.cols = order(cor.mat$cor, decreasing = TRUE),

                        # top plot
                        yt = cor.mat$cor,
                        yt.axis.name = "correlation of\npredicted resp.\nwith true resp.",
                        yt.plot.type = "boxplot",
                        yt.point.size = 4,
                        yt.plot.size = 0.6,
                        yt.axis.size = 15,
                        yt.axis.name.size = 15,
                        #yt.pal = c("#7fcdbb",
                        #           "#41b6c4"),


                        # left labels
                        left.text.angle = 0,
                        left.label.size = 0.25,
                        left.label.text.col = c("black","white"),

                        # bottom labels
                        bottom.label.text.col = c("black","white"),

                        # change to print.plot = TRUE to show the plot
                        print.plot = FALSE)

####### Appendix: fMRI scree plot ####################

# clustered the images
# plotted the responses to the images for the first voxel
# there is a cluster of images who have high responses
voxel <- 1
# combine all voxel responses into a single response (the first prinicpal component)
pc <- princomp(Y)

n <- 10
pc.df <- data.frame(variances = pc$sdev[1:n]^2, component = factor(1:n))
scree <- ggplot(pc.df) +
  geom_bar(aes(x = component, y = variances), stat = "identity") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))


####### Appendix: fMRI number of variables selected ####################

load("fMRI_df.RData")

num_selected <- data.frame(voxel = 1:length(df.lasso), df = df.lasso, cluster = col.membership)
ggplot(dplyr::filter(num_selected, df < 200)) +
  geom_histogram(aes(x = df), col = "white", binwidth = 5) +
  scale_x_continuous(name = "Number of non-zero coefficients") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  facet_wrap(~cluster, ncol = 1) +
  theme(strip.text.x = element_text(size = 16))


num_selected %>% filter(df >= 200) %>% group_by(cluster) %>% summarize(n = n())
num_selected %>% group_by(cluster) %>% summarize(min = min(df), max = max(df))




