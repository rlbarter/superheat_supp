############## cleaning the communities and crimes data ###################


library(igraph)
library(dplyr)

# to_df is a function which transforms an x-y matrix into a data frame
to_df <- function(X) {
  if(!is.matrix(X)) stop("X must be a matrix")

  # converts image matrix into a data frame (one row per observation)
  X.vec <- as.vector(X) # convert the matrix to a vector
  X.mat <- cbind(value = X.vec,
                 # convert vector to matrix with columns
                 # for the x and y coordinates
                 x = rep(1:ncol(X), each = nrow(X)),
                 y = rep(1:nrow(X), times = ncol(X)))
  X.df <- as.data.frame(X.mat)

  return(X.df)
}


# load in the data
crime <- read.delim("../CommViolPredUnnormalizedData.txt", header = F, sep = ",")
names <- as.character(c(read.delim("../CommunityNames.txt", header = F, sep = "\n"))$V1)
colnames(crime) <- names
# remove fold variable
crime <- crime[,-5]


# remove variables with many missing values
to.rm <- which(apply(crime, 2, function(x)  sum(is.na(x))) > 1000)
crime <- crime[,-to.rm]
to.rm <- which(apply(crime, 1, function(x) sum(is.na(x))) > 0)
crime <- crime[-to.rm,]




# remove response varibales
not.var <- c("communityname",
             "state",
             "murders",
             "murdPerPop",
             "rapes",
             "rapesPerPop",
             "robberies",
             "robbbPerPop",
             "assaults",
             "assaultPerPop",
             "burglaries",
             "burglPerPop",
             "larcenies",
             "larcPerPop",
             "autoTheft",
             "autoTheftPerPop",
             "arsons",
             "arsonsPerPop",
             "ViolentCrimesPerPop",
             "nonViolPerPop")

# define the design matrix
X <- crime[,!(colnames(crime) %in% not.var)]
# define the response vector
y <- crime$ViolentCrimesPerPop


################## identify highly correlated clusters of variables ################
# make a graph:
#           nodes are variables
#           edges imply that two variables are highly correlated with one another
cor.df <- to_df(cor(X))
cor.df <- cor.df %>% filter(abs(value) > 0.7)
cor.df <- cor.df %>% filter(!(x == y))
cor.df$x <- colnames(X)[cor.df$x]
cor.df$y <- colnames(X)[cor.df$y]
g <- graph_from_edgelist(as.matrix(cor.df[,-1]), directed = FALSE)
g <- simplify(g)



set.seed(12345)
# decompose into connected subgraphs
g.decom <- decompose(g)
# sample one node/variable from each subgraph
selected.vars <- sapply(g.decom, function(x) sample(as_ids(V(x)),1))
selected.vars <- c("PctKidsBornNeverMar",selected.vars,colnames(X)[!(colnames(X) %in% as_ids(V(g)))])




# restrict to sampled nodes
X.scale <- data.frame(scale(X))
crime.df <- X[,selected.vars]






