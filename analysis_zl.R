### analysis of zl hospital data combining patient and health

setwd("/lustre/user/liclab/shiq/CRC/zl_hospital/analysis")

#############################################################
##### (0)整理并保存patient与health的预处理信息
##########################################################
item.ref.zl = data.frame( item        = patient.ref.zl$item,
                          patient.ref = patient.ref.zl$ref, 
                          health.ref  = health.ref.zl$ref,
                          ref         = patient.ref.zl$ref,
                          stringsAsFactors = FALSE
                          )

item.ref.zl = item.ref.zl[ item.ref.zl$item != "HBeAg", ]
item.ref.zl = item.ref.zl[ item.ref.zl$item != "HBsAg", ]

write.table(item.ref.zl, "/lustre/user/liclab/shiq/CRC/zl_hospital/analysis/item.ref.zl.txt",quote = T)

save.image(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/analysis/zl_hospital_preprocess_result.RData"
  )

#####

#############################################################
##### (1) 预处理加类名
##########################################################
head(patient.info.zl)
patient.matrix.zl <- cbind( patient.info.zl, class = rep( "patient", dim(patient.info.zl)[1] ),
                             stringsAsFactors = FALSE
                             )
head(patient.matrix.zl)
patient.matrix.zl$location[grep("[结]", patient.matrix.zl$location)] <- "colon"
patient.matrix.zl$location[grep("[直]", patient.matrix.zl$location)] <- "rectum"
head(patient.matrix.zl)

head(health.info.zl)
health.matrix.zl <- cbind( health.info.zl, class = rep( "healthy", dim(health.info.zl)[1] ),
                             stringsAsFactors = FALSE
                             )
head(health.matrix.zl)

#####

#############################################################
##### (2) descriptive study of data
##########################################################
################################################ 2.1 table 1
sum( patient.matrix.zl$age<50 ) ## 780
sum( (patient.matrix.zl$age>=50) & (patient.matrix.zl$age<60) ) ## 1164
sum( (patient.matrix.zl$age>=60) & (patient.matrix.zl$age<70) ) ## 1136
sum( (patient.matrix.zl$age>=70) & (patient.matrix.zl$age<80) ) ## 923
sum( patient.matrix.zl$age>=80 ) ## 208
mean( patient.matrix.zl$age ) ## 60.6
sd( patient.matrix.zl$age ) ## 12.5

sum( health.matrix.zl$age<50 ) ## 51220
sum( (health.matrix.zl$age>=50) & (health.matrix.zl$age<60) ) ## 14179
sum( (health.matrix.zl$age>=60) & (health.matrix.zl$age<70) ) ## 6473
sum( (health.matrix.zl$age>=70) & (health.matrix.zl$age<80) ) ## 4319
sum( health.matrix.zl$age>=80 ) ## 908
mean( health.matrix.zl$age ) ## 43.9
sd( health.matrix.zl$age ) ## 14.6

t.test(x=patient.matrix.zl$age, y=health.matrix.zl$age, alternative="greater")
t.test(x=patient.matrix.zl$Alb, y=health.matrix.zl$Alb, alternative="less")
shapiro.test(patient.matrix.zl$age)
wilcox.test(x=patient.matrix.zl$age, y=health.matrix.zl$age)
fisher.test(x=patient.matrix.zl$age, y=health.matrix.zl$age)
anova(x=patient.matrix.zl$age, y=health.matrix.zl$age)
var.test(x=patient.matrix.zl$age, y=health.matrix.zl$age)


table( patient.info.zl$sex )
table( health.info.zl$sex )

table( patient.info.zl$stage )

length( grep("[结]", patient.info.zl$location) ) ## 2279
length( grep("[直]", patient.info.zl$location) ) ## 2057
length( intersect(grep("[结]", patient.info.zl$location), grep("[直]", patient.info.zl$location)) ) ## 125


################################################## 2.2 table 2 p value
pvalue = as.data.frame(matrix(numeric(0),ncol=2))
for(i in colnames( patient.matrix.zl[, -c(1,41)] )){
  pvalue.single = cbind(i, wilcox.test( patient.matrix.zl[,i], health.matrix.zl[,i])$p.value)
  pvalue = rbind( pvalue, pvalue.single)
}
rm(pvalue.single)

################################################## 2.3 fig 1A 1B vioplot
fig.dir = "/lustre/user/liclab/shiq/CRC/zl_hospital/analysis/vioplot/"
library(vioplot)

##### age
png( "age.png", width = 1024, height = 1024 )
par( mar=c(6, 12, 8, 4), cex.axis=4 )

plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=c(0,110),
     axes=FALSE,ann=FALSE
     )
vioplot( patient.matrix.zl$age, health.matrix.zl$age, 
         col = c("orange", "orange"),
         add = T
         )

axis( side=1, at=c(1,2), label=c("CRC", "Healthy"), padj = 1)
axis( side=2, at=seq(0,100,20), las=1 )
mtext( "Age", side=2, line=7, cex=4, padj=0)## las  标签方向
legend("top",
       legend=c("t-test p-value < 2.2e-16"),
       bty="n",
       cex=3)

dev.off()


for(i in colnames(patient.matrix.zl[, -c(1,41)])){

  pvec = patient.matrix.zl[, i][ !is.na( patient.matrix.zl[, i] ) ]
  hvec = health.matrix.zl[, i][ !is.na( health.matrix.zl[, i] ) ]
  
  # unit = units[units$items == i, "units"]
  
#   ref     = item.ref[item.ref.zl$item == i, "ref"] 
#   ref.inf = as.numeric( strsplit( ref, '-')[[1]][1] )
#   ref.sup = as.numeric( strsplit( ref, '-')[[1]][2] )
  
  i = gsub("/", " ratio ", i)
  i = gsub("%", "percent", i)
  png( filename = paste( fig.dir, i, ".png", sep = ''), width = 1024, height = 1024 )
  par( mar=c(8, 4, 4, 4), cex.axis=4 )
  vioplot( pvec, hvec, 
           col = color <- rgb(255,127,0,alpha=80,max=255),
           names = F
  )
  i = gsub(" ratio ", "/", i)
  i = gsub( "percent", "%", i)

  axis( side=1, at=c(1,2), label=c("CRC", "Healthy"), padj = 1)
  # axis( side=2, at=seq(0,100,10), las=1 )
  title( main = list(i, cex=4) )
  # title(  main = i, ylab = unit, cex = 3)
  # abline( h = c( ref.inf, ref.sup ), lty = 2 )
  
  dev.off()
}


############################################ 2.4 sample sizes of every items
item.size.p = apply( patient.matrix.zl[, -41], MARGIN = 2, FUN = function(x){sum(!is.na(x))})
plot(1,1)
barplace = barplot( item.size.p, cex.names = 0.8, srt = 45, axes = FALSE,
                    axisnames = FALSE, ylim = c(-200, 4250)
                    )
axis(2, at = seq(0,4000,500))
text(y = -100, x = barplace, labels =colnames(patient.matrix.zl[, -41]), srt=60, cex = 0.8)

item.size.h = apply( health.matrix.zl[, -41], MARGIN = 2, FUN = function(x){sum(!is.na(x))})
plot(1,1)
barplace = barplot( item.size.h, cex.names = 0.8, srt = 45, axes = FALSE,
                    axisnames = FALSE, ylim = c(-20000, 80000)
)
axis(2, at = c(0, 20000, 40000, 60000, 80000))
text(y = -8000, x = barplace, labels =colnames(health.matrix.zl[, -41]), srt=60, cex = 0.8)

rm(barplace)

#####

#############################################################
##### (3) divide data
##########################################################
##### 3.1 divide all data into test (20%) and train (80%)
dim.zl.p <- dim(patient.matrix.zl)[1]
dim.zl.h <- dim(health.matrix.zl )[1]

dim.zl.p.test <- floor(dim.zl.p*0.2)
dim.zl.h.test <- floor(dim.zl.h*0.2)

set.seed(2000)
data.label.zl.p.test  <- sample( seq(1,dim.zl.p), dim.zl.p.test, replace=FALSE )
data.label.zl.p.train <- setdiff( seq(1,dim.zl.p), data.label.zl.p.test )

set.seed(2000)
data.label.zl.h.test  <- sample( seq(1,dim.zl.h), dim.zl.h.test, replace=FALSE )
data.label.zl.h.train <- setdiff( seq(1,dim.zl.h), data.label.zl.h.test )

data.zl.p.test <- patient.matrix.zl[data.label.zl.p.test, ]
data.zl.h.test <- health.matrix.zl[data.label.zl.h.test, ]

data.zl.p.train <- patient.matrix.zl[data.label.zl.p.train, ]
data.zl.h.train <- health.matrix.zl[data.label.zl.h.train, ]

##### 3.2 divide train data (80%) into valid (20%) and train (60%)
dim.zl.p.train <- dim(data.zl.p.train)[1]
dim.zl.h.train <- dim(data.zl.h.train)[1]

dim.zl.p.valid <- floor(dim.zl.p.train*0.25)
dim.zl.h.valid <- floor(dim.zl.h.train*0.25)

set.seed(3000)
data.label.zl.p.valid <- sample( seq(1,dim.zl.p.train), dim.zl.p.valid, replace=-FALSE )
data.label.zl.p.train <- setdiff( seq(1,dim.zl.p.train), data.label.zl.p.valid )

set.seed(3000)
data.label.zl.h.valid <- sample( seq(1,dim.zl.h.train), dim.zl.h.valid, replace=FALSE )
data.label.zl.h.train <- setdiff( seq(1,dim.zl.h.train), data.label.zl.h.valid )

data.zl.p.valid <- data.zl.p.train[data.label.zl.p.valid, ]
data.zl.h.valid <- data.zl.h.train[data.label.zl.h.valid, ]

data.zl.p.train <- data.zl.p.train[data.label.zl.p.train, ]
data.zl.h.train <- data.zl.h.train[data.label.zl.h.train, ]
data.zl.train   <- rbind( data.zl.p.train[,-c(1,2,3)], data.zl.h.train[,-c(1)] )

#####

#############################################################
##### (4) age model with valid data
##########################################################
tpr.zl.age.valid <- c( sum(data.zl.p.valid$age >= 100)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 95)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 90)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 85)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 80)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 75)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 70)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 65)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 60)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 55)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 50)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 45)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 42)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 40)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 35)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 30)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 25)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 20)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 15)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 10)/length(data.zl.p.valid$age),
                       sum(data.zl.p.valid$age >= 5 )/length(data.zl.p.valid$age)
                       )

fpr.zl.age.valid <- c( sum(data.zl.h.valid$age >= 100)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 95)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 90)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 85)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 80)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 75)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 70)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 65)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 60)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 55)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 50)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 45)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 42)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 40)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 35)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 30)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 25)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 20)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 15)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 10)/length(data.zl.h.valid$age),
                       sum(data.zl.h.valid$age >= 5 )/length(data.zl.h.valid$age)
                       )

#### AUC
AUC.zl.age.valid <- c()
for(i in 2:length(fpr.zl.age.valid)){
  S <-( fpr.zl.age.valid[i]-fpr.zl.age.valid[i-1] ) * 
    ( tpr.zl.age.valid[i-1] + (tpr.zl.age.valid[i]-tpr.zl.age.valid[i-1])/2 )
  AUC.zl.age.valid <- sum(AUC.zl.age.valid, S)
  }

### tpr(fpr=1%)
x1 <- max( fpr.zl.age.valid[fpr.zl.age.valid<=0.01] )
x2 <- min( fpr.zl.age.valid[fpr.zl.age.valid> 0.01] )

y1 <- max( tpr.zl.age.valid[fpr.zl.age.valid==x1] )
y2 <- min( tpr.zl.age.valid[fpr.zl.age.valid==x2] )

tpr.hat.age.valid <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (5) NB training and valid
##########################################################
##################################### 5.1 training
library(e1071)

# fit model
fit.zl.NB <- naiveBayes(class~., data=data.zl.train)

summary(fit.zl.NB)

##################################### 5.2 make predictions
valid.result.zl.p.NB <- predict(fit.zl.NB, data.zl.p.valid[,-c(1,2,3,43)], type="raw")[,2]
valid.result.zl.h.NB <- predict(fit.zl.NB, data.zl.h.valid[,-c(1,41)], type="raw")[,2]

################################### 5.3 计算tpr,fpr
tpr.zl.NB.valid <- c(0)
fpr.zl.NB.valid <- c(0)

valid.result.zl.NB <- unique(c(valid.result.zl.p.NB, valid.result.zl.h.NB))

for(i in sort(valid.result.zl.NB, decreasing=TRUE)){
  tpr.zl.NB.valid <- c(tpr.zl.NB.valid, sum(valid.result.zl.p.NB >= i)/length(valid.result.zl.p.NB))
  fpr.zl.NB.valid <- c(fpr.zl.NB.valid, sum(valid.result.zl.h.NB >= i)/length(valid.result.zl.h.NB))
}

################################## 5.4 计算AUC
AUC.zl.NB.valid <- c()
for(i in 2:length(fpr.zl.NB.valid)){
  S <- ( fpr.zl.NB.valid[i]-fpr.zl.NB.valid[i-1] ) * 
    ( tpr.zl.NB.valid[i-1] + (tpr.zl.NB.valid[i]-tpr.zl.NB.valid[i-1])/2 )
  AUC.zl.NB.valid <- sum(AUC.zl.NB.valid, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.NB.valid[fpr.zl.NB.valid<=0.01] )
x2 <- min( fpr.zl.NB.valid[fpr.zl.NB.valid> 0.01] )

y1 <- max( tpr.zl.NB.valid[fpr.zl.NB.valid==x1] )
y2 <- min( tpr.zl.NB.valid[fpr.zl.NB.valid==x2] )

tpr.hat.NB.valid <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (6) CART training and valid
##########################################################
library(rpart)

### select best CART

##################################### 6.1 training
fit.zl.CART <- rpart(class~., data=data.zl.train[,c("age", "Alb", "HCT", "LYMPH%", "class")])

summary(fit.zl.CART)

##################################### 6.2 make predictions
# make predictions
valid.result.zl.p.CART <- predict(fit.zl.CART, data.zl.p.valid[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
valid.result.zl.h.CART <- predict(fit.zl.CART, data.zl.h.valid[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]

################################### 6.3 计算tpr,fpr
tpr.zl.CART.valid <- c(0)
fpr.zl.CART.valid <- c(0)

valid.result.zl.CART <- unique(c(valid.result.zl.p.CART, valid.result.zl.h.CART))

for(i in sort(valid.result.zl.CART, decreasing=TRUE)){
  tpr.zl.CART.valid <- c(tpr.zl.CART.valid, sum(valid.result.zl.p.CART >= i)/length(valid.result.zl.p.CART))
  fpr.zl.CART.valid <- c(fpr.zl.CART.valid, sum(valid.result.zl.h.CART >= i)/length(valid.result.zl.h.CART))
}

################################## 5.4 计算AUC
AUC.zl.CART.valid <- c()
for(i in 2:length(fpr.zl.CART.valid)){
  S <- ( fpr.zl.CART.valid[i]-fpr.zl.CART.valid[i-1] ) * 
    ( tpr.zl.CART.valid[i-1] + (tpr.zl.CART.valid[i]-tpr.zl.CART.valid[i-1])/2 )
  AUC.zl.CART.valid <- sum(AUC.zl.CART.valid, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.CART.valid[fpr.zl.CART.valid<=0.01] )
x2 <- min( fpr.zl.CART.valid[fpr.zl.CART.valid> 0.01] )

y1 <- max( tpr.zl.CART.valid[fpr.zl.CART.valid==x1] )
y2 <- min( tpr.zl.CART.valid[fpr.zl.CART.valid==x2] )

tpr.hat.CART.valid <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (7) age model with testing data
##########################################################
tpr.zl.age.test <- c( sum(data.zl.p.test$age >= 100)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 95)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 90)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 85)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 80)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 75)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 70)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 65)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 60)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 55)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 50)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 45)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 42)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 40)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 35)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 30)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 25)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 20)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 15)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 10)/length(data.zl.p.test$age),
                      sum(data.zl.p.test$age >= 5 )/length(data.zl.p.test$age)
                      )

fpr.zl.age.test <- c( sum(data.zl.h.test$age >= 100)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 95)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 90)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 85)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 80)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 75)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 70)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 65)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 60)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 55)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 50)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 45)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 42)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 40)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 35)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 30)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 25)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 20)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 15)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 10)/length(data.zl.h.test$age),
                      sum(data.zl.h.test$age >= 5 )/length(data.zl.h.test$age)
                      )

### AUC
AUC.zl.age.test <- c()
for(i in 2:length(fpr.zl.age.test)){
  S <- ( fpr.zl.age.test[i]-fpr.zl.age.test[i-1] ) * 
    ( tpr.zl.age.test[i-1] + (tpr.zl.age.test[i]-tpr.zl.age.test[i-1])/2 )
  AUC.zl.age.test <- sum(AUC.zl.age.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.age.test[fpr.zl.age.test<=0.01] )
x2 <- min( fpr.zl.age.test[fpr.zl.age.test> 0.01] )

y1 <- max( tpr.zl.age.test[fpr.zl.age.test==x1] )
y2 <- min( tpr.zl.age.test[fpr.zl.age.test==x2] )

tpr.hat.age.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (8) NB testing
##########################################################
##################################### 8.1 predicting
test.result.zl.p.NB <- predict(fit.zl.NB, data.zl.p.test[,-c(1,2,3,43)], type="raw")[,2]
test.result.zl.h.NB <- predict(fit.zl.NB, data.zl.h.test[,-c(1,41)], type="raw")[,2]

################################### 8.2 计算tpr,fpr
tpr.zl.NB.test <- c(0)
fpr.zl.NB.test <- c(0)

test.result.zl.NB <- unique(c(test.result.zl.p.NB, test.result.zl.h.NB))

for(i in sort(test.result.zl.NB, decreasing=TRUE)){
  tpr.zl.NB.test <- c(tpr.zl.NB.test, sum(test.result.zl.p.NB >= i)/length(test.result.zl.p.NB))
  fpr.zl.NB.test <- c(fpr.zl.NB.test, sum(test.result.zl.h.NB >= i)/length(test.result.zl.h.NB))
}

################################## 8.3 计算AUC
AUC.zl.NB.test <- c()
for(i in 2:length(fpr.zl.NB.test)){
  S <- ( fpr.zl.NB.test[i]-fpr.zl.NB.test[i-1] ) * 
    ( tpr.zl.NB.test[i-1] + (tpr.zl.NB.test[i]-tpr.zl.NB.test[i-1])/2 )
  AUC.zl.NB.test <- sum(AUC.zl.NB.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.NB.test[fpr.zl.NB.test<=0.01] )
x2 <- min( fpr.zl.NB.test[fpr.zl.NB.test> 0.01] )

y1 <- max( tpr.zl.NB.test[fpr.zl.NB.test==x1] )
y2 <- min( tpr.zl.NB.test[fpr.zl.NB.test==x2] )

tpr.hat.NB.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (9) CART testing
##########################################################
##################################### 9.1 predicting
test.result.zl.p.CART <- predict(fit.zl.CART, data.zl.p.test[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.zl.h.CART <- predict(fit.zl.CART, data.zl.h.test[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]

################################### 9.2 计算tpr,fpr
tpr.zl.CART.test <- c(0)
fpr.zl.CART.test <- c(0)

test.result.zl.CART <- unique(c(test.result.zl.p.CART, test.result.zl.h.CART))

for(i in sort(test.result.zl.CART, decreasing=TRUE)){
  tpr.zl.CART.test <- c(tpr.zl.CART.test, sum(test.result.zl.p.CART >= i)/length(test.result.zl.p.CART))
  fpr.zl.CART.test <- c(fpr.zl.CART.test, sum(test.result.zl.h.CART >= i)/length(test.result.zl.h.CART))
}

################################## 9.3 计算AUC
AUC.zl.CART.test <- c()
for(i in 2:length(fpr.zl.CART.test)){
  S <- ( fpr.zl.CART.test[i]-fpr.zl.CART.test[i-1] ) * 
    ( tpr.zl.CART.test[i-1] + (tpr.zl.CART.test[i]-tpr.zl.CART.test[i-1])/2 )
  AUC.zl.CART.test <- sum(AUC.zl.CART.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.CART.test[fpr.zl.CART.test<=0.01] )
x2 <- min( fpr.zl.CART.test[fpr.zl.CART.test> 0.01] )

y1 <- max( tpr.zl.CART.test[fpr.zl.CART.test==x1] )
y2 <- min( tpr.zl.CART.test[fpr.zl.CART.test==x2] )

tpr.hat.CART.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (10) CART test (ratio)
##########################################################
#################################### rate01 800 : 10000
AUC.zl.test.ratio01 <- c()
tpr.hat.ratio01 <- c()

for(j in 1:100){
###### 1 组合样本
set.seed(j)
test.result.zl.p.ratio01 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),800,replace=FALSE)]
set.seed(j)
test.result.zl.h.ratio01 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),10000,replace=FALSE)]

####### 2 fpr,tpr
tpr.zl.test.ratio01 <- c(0)
fpr.zl.test.ratio01 <- c(0)

test.result.zl.ratio01 <- unique(c(test.result.zl.p.ratio01, test.result.zl.h.ratio01))

for(i in sort(test.result.zl.ratio01, decreasing=TRUE)){
  tpr.zl.test.ratio01 <- c(tpr.zl.test.ratio01, 
                          sum(test.result.zl.p.ratio01 >= i)/length(test.result.zl.p.ratio01))
  fpr.zl.test.ratio01 <- c(fpr.zl.test.ratio01,
                          sum(test.result.zl.h.ratio01 >= i)/length(test.result.zl.h.ratio01))
}

####### 3 计算AUC
AUC.temp <- c()
for(i in 2:length(fpr.zl.test.ratio01)){
  S = ( fpr.zl.test.ratio01[i] - fpr.zl.test.ratio01[i-1] ) * 
    ( tpr.zl.test.ratio01[i-1] + (tpr.zl.test.ratio01[i] - tpr.zl.test.ratio01[i-1])/2 )
  AUC.temp <- sum(AUC.temp, S)
}
AUC.zl.test.ratio01 = c(AUC.zl.test.ratio01, AUC.temp)

####### 4 计算tpr(fpr=1%)
x1 <- max(fpr.zl.test.ratio01[fpr.zl.test.ratio01<=0.01])
x2 <- min(fpr.zl.test.ratio01[fpr.zl.test.ratio01> 0.01])

y1 <- tpr.zl.test.ratio01[fpr.zl.test.ratio01==x1]
y2 <- tpr.zl.test.ratio01[fpr.zl.test.ratio01==x2]

tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
tpr.hat.ratio01 <- c(tpr.hat.ratio01, tpr.hat.temp)

}

se.AUC.zl.test.ratio01 <- sd((AUC.zl.test.ratio01))/10
AUC.zl.test.ratio01 <- mean(AUC.zl.test.ratio01)

se.tpr.hat.ratio01 <- sd(tpr.hat.ratio01)/10
tpr.hat.ratio01 <- mean(tpr.hat.ratio01)

#################################### rate02 800 : 15000
AUC.zl.test.ratio02 <- c()
tpr.hat.ratio02 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(j)
  test.result.zl.p.ratio02 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),800,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio02 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio02 <- c(0)
  fpr.zl.test.ratio02 <- c(0)
  
  test.result.zl.ratio02 <- unique(c(test.result.zl.p.ratio02, test.result.zl.h.ratio02))
  
  for(i in sort(test.result.zl.ratio02, decreasing=TRUE)){
    tpr.zl.test.ratio02 <- c(tpr.zl.test.ratio02, 
                             sum(test.result.zl.p.ratio02 >= i)/length(test.result.zl.p.ratio02))
    fpr.zl.test.ratio02 <- c(fpr.zl.test.ratio02,
                             sum(test.result.zl.h.ratio02 >= i)/length(test.result.zl.h.ratio02))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio02)){
    S = ( fpr.zl.test.ratio02[i] - fpr.zl.test.ratio02[i-1] ) * 
      ( tpr.zl.test.ratio02[i-1] + (tpr.zl.test.ratio02[i] - tpr.zl.test.ratio02[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio02 = c(AUC.zl.test.ratio02, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio02[fpr.zl.test.ratio02<=0.01])
  x2 <- min(fpr.zl.test.ratio02[fpr.zl.test.ratio02> 0.01])
  
  y1 <- tpr.zl.test.ratio02[fpr.zl.test.ratio02==x1]
  y2 <- tpr.zl.test.ratio02[fpr.zl.test.ratio02==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio02 <- c(tpr.hat.ratio02, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio02 <- sd((AUC.zl.test.ratio02))/10
AUC.zl.test.ratio02 <- mean(AUC.zl.test.ratio02)

se.tpr.hat.ratio02 <- sd(tpr.hat.ratio02)/10
tpr.hat.ratio02 <- mean(tpr.hat.ratio02)

#################################### rate03 500 : 15000
AUC.zl.test.ratio03 <- c()
tpr.hat.ratio03 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(j)
  test.result.zl.p.ratio03 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),500,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio03 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio03 <- c(0)
  fpr.zl.test.ratio03 <- c(0)
  
  test.result.zl.ratio03 <- unique(c(test.result.zl.p.ratio03, test.result.zl.h.ratio03))
  
  for(i in sort(test.result.zl.ratio03, decreasing=TRUE)){
    tpr.zl.test.ratio03 <- c(tpr.zl.test.ratio03, 
                             sum(test.result.zl.p.ratio03 >= i)/length(test.result.zl.p.ratio03))
    fpr.zl.test.ratio03 <- c(fpr.zl.test.ratio03,
                             sum(test.result.zl.h.ratio03 >= i)/length(test.result.zl.h.ratio03))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio03)){
    S = ( fpr.zl.test.ratio03[i] - fpr.zl.test.ratio03[i-1] ) * 
      ( tpr.zl.test.ratio03[i-1] + (tpr.zl.test.ratio03[i] - tpr.zl.test.ratio03[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio03 = c(AUC.zl.test.ratio03, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio03[fpr.zl.test.ratio03<=0.01])
  x2 <- min(fpr.zl.test.ratio03[fpr.zl.test.ratio03> 0.01])
  
  y1 <- tpr.zl.test.ratio03[fpr.zl.test.ratio03==x1]
  y2 <- tpr.zl.test.ratio03[fpr.zl.test.ratio03==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio03 <- c(tpr.hat.ratio03, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio03 <- sd((AUC.zl.test.ratio03))/10
AUC.zl.test.ratio03 <- mean(AUC.zl.test.ratio03)

se.tpr.hat.ratio03 <- sd(tpr.hat.ratio03)/10
tpr.hat.ratio03 <- mean(tpr.hat.ratio03)

#################################### rate04 300 : 15000
AUC.zl.test.ratio04 <- c()
tpr.hat.ratio04 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(j)
  test.result.zl.p.ratio04 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),300,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio04 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio04 <- c(0)
  fpr.zl.test.ratio04 <- c(0)
  
  test.result.zl.ratio04 <- unique(c(test.result.zl.p.ratio04, test.result.zl.h.ratio04))
  
  for(i in sort(test.result.zl.ratio04, decreasing=TRUE)){
    tpr.zl.test.ratio04 <- c(tpr.zl.test.ratio04, 
                             sum(test.result.zl.p.ratio04 >= i)/length(test.result.zl.p.ratio04))
    fpr.zl.test.ratio04 <- c(fpr.zl.test.ratio04,
                             sum(test.result.zl.h.ratio04 >= i)/length(test.result.zl.h.ratio04))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio04)){
    S = ( fpr.zl.test.ratio04[i] - fpr.zl.test.ratio04[i-1] ) * 
      ( tpr.zl.test.ratio04[i-1] + (tpr.zl.test.ratio04[i] - tpr.zl.test.ratio04[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio04 = c(AUC.zl.test.ratio04, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio04[fpr.zl.test.ratio04<=0.01])
  x2 <- min(fpr.zl.test.ratio04[fpr.zl.test.ratio04> 0.01])
  
  y1 <- tpr.zl.test.ratio04[fpr.zl.test.ratio04==x1]
  y2 <- tpr.zl.test.ratio04[fpr.zl.test.ratio04==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio04 <- c(tpr.hat.ratio04, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio04 <- sd((AUC.zl.test.ratio04))/10
AUC.zl.test.ratio04 <- mean(AUC.zl.test.ratio04)

se.tpr.hat.ratio04 <- sd(tpr.hat.ratio04)/10
tpr.hat.ratio04 <- mean(tpr.hat.ratio04)

#################################### rate05 100 : 15000
AUC.zl.test.ratio05 <- c()
tpr.hat.ratio05 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(j)
  test.result.zl.p.ratio05 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),100,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio05 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio05 <- c(0)
  fpr.zl.test.ratio05 <- c(0)
  
  test.result.zl.ratio05 <- unique(c(test.result.zl.p.ratio05, test.result.zl.h.ratio05))
  
  for(i in sort(test.result.zl.ratio05, decreasing=TRUE)){
    tpr.zl.test.ratio05 <- c(tpr.zl.test.ratio05, 
                             sum(test.result.zl.p.ratio05 >= i)/length(test.result.zl.p.ratio05))
    fpr.zl.test.ratio05 <- c(fpr.zl.test.ratio05,
                             sum(test.result.zl.h.ratio05 >= i)/length(test.result.zl.h.ratio05))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio05)){
    S = ( fpr.zl.test.ratio05[i] - fpr.zl.test.ratio05[i-1] ) * 
      ( tpr.zl.test.ratio05[i-1] + (tpr.zl.test.ratio05[i] - tpr.zl.test.ratio05[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio05 = c(AUC.zl.test.ratio05, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio05[fpr.zl.test.ratio05<=0.01])
  x2 <- min(fpr.zl.test.ratio05[fpr.zl.test.ratio05> 0.01])
  
  y1 <- tpr.zl.test.ratio05[fpr.zl.test.ratio05==x1]
  y2 <- tpr.zl.test.ratio05[fpr.zl.test.ratio05==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio05 <- c(tpr.hat.ratio05, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio05 <- sd((AUC.zl.test.ratio05))/10
AUC.zl.test.ratio05 <- mean(AUC.zl.test.ratio05)

se.tpr.hat.ratio05 <- sd(tpr.hat.ratio05)/10
tpr.hat.ratio05 <- mean(tpr.hat.ratio05)

#################################### rate06 50 : 15000
AUC.zl.test.ratio06 <- c()
tpr.hat.ratio06 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(j)
  test.result.zl.p.ratio06 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),50,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio06 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio06 <- c(0)
  fpr.zl.test.ratio06 <- c(0)
  
  test.result.zl.ratio06 <- unique(c(test.result.zl.p.ratio06, test.result.zl.h.ratio06))
  
  for(i in sort(test.result.zl.ratio06, decreasing=TRUE)){
    tpr.zl.test.ratio06 <- c(tpr.zl.test.ratio06, 
                             sum(test.result.zl.p.ratio06 >= i)/length(test.result.zl.p.ratio06))
    fpr.zl.test.ratio06 <- c(fpr.zl.test.ratio06,
                             sum(test.result.zl.h.ratio06 >= i)/length(test.result.zl.h.ratio06))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio06)){
    S = ( fpr.zl.test.ratio06[i] - fpr.zl.test.ratio06[i-1] ) * 
      ( tpr.zl.test.ratio06[i-1] + (tpr.zl.test.ratio06[i] - tpr.zl.test.ratio06[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio06 = c(AUC.zl.test.ratio06, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio06[fpr.zl.test.ratio06<=0.01])
  x2 <- min(fpr.zl.test.ratio06[fpr.zl.test.ratio06> 0.01])
  
  y1 <- tpr.zl.test.ratio06[fpr.zl.test.ratio06==x1]
  y2 <- tpr.zl.test.ratio06[fpr.zl.test.ratio06==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio06 <- c(tpr.hat.ratio06, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio06 <- sd((AUC.zl.test.ratio06))/10
AUC.zl.test.ratio06 <- mean(AUC.zl.test.ratio06)

se.tpr.hat.ratio06 <- sd(tpr.hat.ratio06)/10
tpr.hat.ratio06 <- mean(tpr.hat.ratio06)


#################################### rate07 30 : 15000
AUC.zl.test.ratio07 <- c()
tpr.hat.ratio07 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(j)
  test.result.zl.p.ratio07 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),30,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio07 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio07 <- c(0)
  fpr.zl.test.ratio07 <- c(0)
  
  test.result.zl.ratio07 <- unique(c(test.result.zl.p.ratio07, test.result.zl.h.ratio07))
  
  for(i in sort(test.result.zl.ratio07, decreasing=TRUE)){
    tpr.zl.test.ratio07 <- c(tpr.zl.test.ratio07, 
                             sum(test.result.zl.p.ratio07 >= i)/length(test.result.zl.p.ratio07))
    fpr.zl.test.ratio07 <- c(fpr.zl.test.ratio07,
                             sum(test.result.zl.h.ratio07 >= i)/length(test.result.zl.h.ratio07))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio07)){
    S = ( fpr.zl.test.ratio07[i] - fpr.zl.test.ratio07[i-1] ) * 
      ( tpr.zl.test.ratio07[i-1] + (tpr.zl.test.ratio07[i] - tpr.zl.test.ratio07[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio07 = c(AUC.zl.test.ratio07, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio07[fpr.zl.test.ratio07<=0.01])
  x2 <- min(fpr.zl.test.ratio07[fpr.zl.test.ratio07> 0.01])
  
  y1 <- tpr.zl.test.ratio07[fpr.zl.test.ratio07==x1]
  y2 <- tpr.zl.test.ratio07[fpr.zl.test.ratio07==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio07 <- c(tpr.hat.ratio07, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio07 <- sd((AUC.zl.test.ratio07))/10
AUC.zl.test.ratio07 <- mean(AUC.zl.test.ratio07)

se.tpr.hat.ratio07 <- sd(tpr.hat.ratio07)/10
tpr.hat.ratio07 <- mean(tpr.hat.ratio07)

#################################### rate08 10 : 15000
AUC.zl.test.ratio08 <- c()
tpr.hat.ratio08 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(j)
  test.result.zl.p.ratio08 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),10,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio08 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio08 <- c(0)
  fpr.zl.test.ratio08 <- c(0)
  
  test.result.zl.ratio08 <- unique(c(test.result.zl.p.ratio08, test.result.zl.h.ratio08))
  
  for(i in sort(test.result.zl.ratio08, decreasing=TRUE)){
    tpr.zl.test.ratio08 <- c(tpr.zl.test.ratio08, 
                             sum(test.result.zl.p.ratio08 >= i)/length(test.result.zl.p.ratio08))
    fpr.zl.test.ratio08 <- c(fpr.zl.test.ratio08,
                             sum(test.result.zl.h.ratio08 >= i)/length(test.result.zl.h.ratio08))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio08)){
    S = ( fpr.zl.test.ratio08[i] - fpr.zl.test.ratio08[i-1] ) * 
      ( tpr.zl.test.ratio08[i-1] + (tpr.zl.test.ratio08[i] - tpr.zl.test.ratio08[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio08 = c(AUC.zl.test.ratio08, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio08[fpr.zl.test.ratio08<=0.01])
  x2 <- min(fpr.zl.test.ratio08[fpr.zl.test.ratio08> 0.01])
  
  y1 <- tpr.zl.test.ratio08[fpr.zl.test.ratio08==x1]
  y2 <- tpr.zl.test.ratio08[fpr.zl.test.ratio08==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio08 <- c(tpr.hat.ratio08, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio08 <- sd((AUC.zl.test.ratio08))/10
AUC.zl.test.ratio08 <- mean(AUC.zl.test.ratio08)

se.tpr.hat.ratio08 <- sd(tpr.hat.ratio08)/10
tpr.hat.ratio08 <- mean(tpr.hat.ratio08)

#################################### rate09 4 : 15000
AUC.zl.test.ratio09 <- c()
tpr.hat.ratio09 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(j)
  test.result.zl.p.ratio09 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),4,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio09 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio09 <- c(0)
  fpr.zl.test.ratio09 <- c(0)
  
  test.result.zl.ratio09 <- unique(c(test.result.zl.p.ratio09, test.result.zl.h.ratio09))
  
  for(i in sort(test.result.zl.ratio09, decreasing=TRUE)){
    tpr.zl.test.ratio09 <- c(tpr.zl.test.ratio09, 
                             sum(test.result.zl.p.ratio09 >= i)/length(test.result.zl.p.ratio09))
    fpr.zl.test.ratio09 <- c(fpr.zl.test.ratio09,
                             sum(test.result.zl.h.ratio09 >= i)/length(test.result.zl.h.ratio09))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio09)){
    S = ( fpr.zl.test.ratio09[i] - fpr.zl.test.ratio09[i-1] ) * 
      ( tpr.zl.test.ratio09[i-1] + (tpr.zl.test.ratio09[i] - tpr.zl.test.ratio09[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio09 = c(AUC.zl.test.ratio09, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio09[fpr.zl.test.ratio09<=0.01])
  x2 <- min(fpr.zl.test.ratio09[fpr.zl.test.ratio09> 0.01])
  
  y1 <- tpr.zl.test.ratio09[fpr.zl.test.ratio09==x1]
  y2 <- tpr.zl.test.ratio09[fpr.zl.test.ratio09==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio09 <- c(tpr.hat.ratio09, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio09 <- sd((AUC.zl.test.ratio09))/10
AUC.zl.test.ratio09 <- mean(AUC.zl.test.ratio09)

se.tpr.hat.ratio09 <- sd(tpr.hat.ratio09)/10
tpr.hat.ratio09 <- mean(tpr.hat.ratio09)

#################################### rate10 2 : 15000
AUC.zl.test.ratio10 <- c()
tpr.hat.ratio10 <- c()

for(j in 1:100){
  ###### 1 组合样本
  set.seed(2*j)
  test.result.zl.p.ratio10 <- test.result.zl.p.CART[sample(1:length(test.result.zl.p.CART),2,replace=FALSE)]
  set.seed(j)
  test.result.zl.h.ratio10 <- test.result.zl.h.CART[sample(1:length(test.result.zl.h.CART),15000,replace=FALSE)]
  
  ####### 2 fpr,tpr
  tpr.zl.test.ratio10 <- c(0)
  fpr.zl.test.ratio10 <- c(0)
  
  test.result.zl.ratio10 <- unique(c(test.result.zl.p.ratio10, test.result.zl.h.ratio10))
  
  for(i in sort(test.result.zl.ratio10, decreasing=TRUE)){
    tpr.zl.test.ratio10 <- c(tpr.zl.test.ratio10, 
                             sum(test.result.zl.p.ratio10 >= i)/length(test.result.zl.p.ratio10))
    fpr.zl.test.ratio10 <- c(fpr.zl.test.ratio10,
                             sum(test.result.zl.h.ratio10 >= i)/length(test.result.zl.h.ratio10))
  }
  
  ####### 3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.test.ratio10)){
    S = ( fpr.zl.test.ratio10[i] - fpr.zl.test.ratio10[i-1] ) * 
      ( tpr.zl.test.ratio10[i-1] + (tpr.zl.test.ratio10[i] - tpr.zl.test.ratio10[i-1])/2 )
    AUC.temp <- sum(AUC.temp, S)
  }
  AUC.zl.test.ratio10 = c(AUC.zl.test.ratio10, AUC.temp)
  
  ####### 4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.test.ratio10[fpr.zl.test.ratio10<=0.01])
  x2 <- min(fpr.zl.test.ratio10[fpr.zl.test.ratio10> 0.01])
  
  y1 <- tpr.zl.test.ratio10[fpr.zl.test.ratio10==x1]
  y2 <- tpr.zl.test.ratio10[fpr.zl.test.ratio10==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.ratio10 <- c(tpr.hat.ratio10, tpr.hat.temp)
  
}

se.AUC.zl.test.ratio10 <- sd((AUC.zl.test.ratio10))/10
AUC.zl.test.ratio10 <- mean(AUC.zl.test.ratio10)

se.tpr.hat.ratio10 <- sd(tpr.hat.ratio10)/10
tpr.hat.ratio10 <- mean(tpr.hat.ratio10)



#### 处理

AUC.zl.test.ratio <- c(AUC.zl.test.ratio01, AUC.zl.test.ratio02, AUC.zl.test.ratio03,
                       AUC.zl.test.ratio04, AUC.zl.test.ratio05, AUC.zl.test.ratio06,
                       AUC.zl.test.ratio07, AUC.zl.test.ratio08, AUC.zl.test.ratio09,
                       AUC.zl.test.ratio10
                       )

se.AUC.zl.test.ratio <- c(se.AUC.zl.test.ratio01, se.AUC.zl.test.ratio02, se.AUC.zl.test.ratio03,
                          se.AUC.zl.test.ratio04, se.AUC.zl.test.ratio05, se.AUC.zl.test.ratio06,
                          se.AUC.zl.test.ratio07, se.AUC.zl.test.ratio08, se.AUC.zl.test.ratio09,
                          se.AUC.zl.test.ratio10
                          )

errorbar.ratio.AUC.low  <- AUC.zl.test.ratio - se.AUC.zl.test.ratio
errorbar.ratio.AUC.high <- AUC.zl.test.ratio + se.AUC.zl.test.ratio


tpr.hat.ratio <- c(tpr.hat.ratio01, tpr.hat.ratio02, tpr.hat.ratio03,
                   tpr.hat.ratio04, tpr.hat.ratio05, tpr.hat.ratio06,
                   tpr.hat.ratio07, tpr.hat.ratio08, tpr.hat.ratio09,
                   tpr.hat.ratio10
                   )

se.tpr.hat.ratio <- c(se.tpr.hat.ratio01, se.tpr.hat.ratio02, se.tpr.hat.ratio03,
                      se.tpr.hat.ratio04, se.tpr.hat.ratio05, se.tpr.hat.ratio06,
                      se.tpr.hat.ratio07, se.tpr.hat.ratio08, se.tpr.hat.ratio09,
                      se.tpr.hat.ratio09
                      )

errorbar.ratio.tpr.low  <- tpr.hat.ratio - se.tpr.hat.ratio
errorbar.ratio.tpr.high <- tpr.hat.ratio + se.tpr.hat.ratio



######

#############################################################
##### (11) CART test (老龄化)
##########################################################
################################ 10.1  分组age>60
####      >60      <=60
#### P    420       422
#### H    2079     13340
################################ 10.1  分组age>60

data.zl.p.test.60o  <- data.zl.p.test[data.zl.p.test$age > 60, c("age", "Alb", "HCT", "LYMPH%", "class")]
data.zl.p.test.60y  <- data.zl.p.test[data.zl.p.test$age <=60, c("age", "Alb", "HCT", "LYMPH%", "class")]

data.zl.h.test.60o  <- data.zl.h.test[data.zl.h.test$age > 60, c("age", "Alb", "HCT", "LYMPH%", "class")]
data.zl.h.test.60y  <- data.zl.h.test[data.zl.h.test$age <=60, c("age", "Alb", "HCT", "LYMPH%", "class")]


################################ 10.2 make predictions
test.result.zl.p.60o <- predict(fit.zl.CART, data.zl.p.test.60o[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.zl.p.60y <- predict(fit.zl.CART, data.zl.p.test.60y[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]

test.result.zl.h.60o <- predict(fit.zl.CART, data.zl.h.test.60o[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.zl.h.60y <- predict(fit.zl.CART, data.zl.h.test.60y[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]


################################### 10.3 0%
AUC.zl.60.test.0perc <- c()
tpr.hat.0perc <- c()

for(j in 1:100){
####### 10.3.1 组合样本
test.result.zl.60.p.0perc <- test.result.zl.p.60y
test.result.zl.60.h.0perc <- test.result.zl.h.60y

####### 10.3.2 fpr,tpr
tpr.zl.60.test.0perc <- c(0)
fpr.zl.60.test.0perc <- c(0)

test.result.zl.60.0perc = unique(c(test.result.zl.60.p.0perc, test.result.zl.60.h.0perc))

for(i in sort(test.result.zl.60.0perc, decreasing = T)){
  tpr.zl.60.test.0perc = c(tpr.zl.60.test.0perc, 
                           sum(test.result.zl.60.p.0perc >= i)/length(test.result.zl.60.p.0perc))
  fpr.zl.60.test.0perc = c(fpr.zl.60.test.0perc,
                           sum(test.result.zl.60.h.0perc >= i)/length(test.result.zl.60.h.0perc))
  }

####### 10.3.3 计算AUC
AUC.temp <- c()
for(i in 2:length(fpr.zl.60.test.0perc)){
  S = ( fpr.zl.60.test.0perc[i] - fpr.zl.60.test.0perc[i-1] ) * 
    ( tpr.zl.60.test.0perc[i-1] + (tpr.zl.60.test.0perc[i] - tpr.zl.60.test.0perc[i-1])/2 )
  AUC.temp <- sum(AUC.temp, S)
}
AUC.zl.60.test.0perc <- c(AUC.zl.60.test.0perc, AUC.temp)

####### 10.3.4 计算tpr(fpr=1%)
x1 <- max(fpr.zl.60.test.0perc[fpr.zl.60.test.0perc<=0.01])
x2 <- min(fpr.zl.60.test.0perc[fpr.zl.60.test.0perc> 0.01])

y1 <- tpr.zl.60.test.0perc[fpr.zl.60.test.0perc==x1]
y2 <- tpr.zl.60.test.0perc[fpr.zl.60.test.0perc==x2]

tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
tpr.hat.0perc <- c(tpr.hat.0perc, tpr.hat.temp)

}

se.AUC.zl.60.test.0perc <- sd(AUC.zl.60.test.0perc)/10
AUC.zl.60.test.0perc <- mean(AUC.zl.60.test.0perc)

se.tpr.hat.0perc <- sd(tpr.hat.0perc)/10
tpr.hat.0perc <- mean(tpr.hat.0perc)


################################### 10.4 5%
AUC.zl.60.test.5perc <- c()
tpr.hat.5perc <- c()

for(j in 1:100){
####### 10.4.1 组合样本
set.seed(j)
temp.60o.p <- test.result.zl.p.60o[sample(1:length(test.result.zl.p.60o),22,replace=FALSE)]
test.result.zl.60.p.5perc <- c(test.result.zl.p.60y, temp.60o.p)

set.seed(j)
temp.60o.h <- test.result.zl.h.60o[sample(1:length(test.result.zl.h.60o),702,replace=FALSE)]
test.result.zl.60.h.5perc <- c(test.result.zl.h.60y, temp.60o.h)

####### 10.4.2 fpr,tpr
tpr.zl.60.test.5perc <- c(0)
fpr.zl.60.test.5perc <- c(0)

test.result.zl.60.5perc = unique(c(test.result.zl.60.p.5perc, test.result.zl.60.h.5perc))

for(i in sort(test.result.zl.60.5perc, decreasing = T)){
  tpr.zl.60.test.5perc = c(tpr.zl.60.test.5perc, 
                           sum(test.result.zl.60.p.5perc >= i)/length(test.result.zl.60.p.5perc))
  fpr.zl.60.test.5perc = c(fpr.zl.60.test.5perc,
                           sum(test.result.zl.60.h.5perc >= i)/length(test.result.zl.60.h.5perc))
}

####### 10.4.3 计算AUC
AUC.temp <- c()
for(i in 2:length(fpr.zl.60.test.5perc)){
  S = ( fpr.zl.60.test.5perc[i] - fpr.zl.60.test.5perc[i-1] ) * 
    ( tpr.zl.60.test.5perc[i-1] + (tpr.zl.60.test.5perc[i] - tpr.zl.60.test.5perc[i-1])/2 )
  AUC.temp = sum(AUC.temp, S)
}
AUC.zl.60.test.5perc <- c(AUC.zl.60.test.5perc, AUC.temp)

####### 10.4.4 计算tpr(fpr=1%)
x1 <- max(fpr.zl.60.test.5perc[fpr.zl.60.test.5perc<=0.01])
x2 <- min(fpr.zl.60.test.5perc[fpr.zl.60.test.5perc> 0.01])

y1 <- tpr.zl.60.test.5perc[fpr.zl.60.test.5perc==x1]
y2 <- tpr.zl.60.test.5perc[fpr.zl.60.test.5perc==x2]

tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
tpr.hat.5perc <- c(tpr.hat.5perc, tpr.hat.temp)

}

se.AUC.zl.60.test.5perc <- sd(AUC.zl.60.test.5perc)/10
AUC.zl.60.test.5perc <- mean(AUC.zl.60.test.5perc)

se.tpr.hat.5perc <- sd(tpr.hat.5perc)/10
tpr.hat.5perc <- mean(tpr.hat.5perc)


################################### 10.5 10%
AUC.zl.60.test.10perc <- c()
tpr.hat.10perc <- c()

for(j in 1:100){
  ####### 10.4.1 组合样本
  set.seed(j)
  temp.60o.p <- test.result.zl.p.60o[sample(1:length(test.result.zl.p.60o),47,replace=FALSE)]
  test.result.zl.60.p.10perc <- c(test.result.zl.p.60y, temp.60o.p)
  
  set.seed(j)
  temp.60o.h <- test.result.zl.h.60o[sample(1:length(test.result.zl.h.60o),1482,replace=FALSE)]
  test.result.zl.60.h.10perc <- c(test.result.zl.h.60y, temp.60o.h)
  
  ####### 10.4.2 fpr,tpr
  tpr.zl.60.test.10perc <- c(0)
  fpr.zl.60.test.10perc <- c(0)
  
  test.result.zl.60.10perc = unique(c(test.result.zl.60.p.10perc, test.result.zl.60.h.10perc))
  
  for(i in sort(test.result.zl.60.10perc, decreasing = T)){
    tpr.zl.60.test.10perc = c(tpr.zl.60.test.10perc, 
                             sum(test.result.zl.60.p.10perc >= i)/length(test.result.zl.60.p.10perc))
    fpr.zl.60.test.10perc = c(fpr.zl.60.test.10perc,
                             sum(test.result.zl.60.h.10perc >= i)/length(test.result.zl.60.h.10perc))
  }
  
  ####### 10.4.3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.60.test.10perc)){
    S = ( fpr.zl.60.test.10perc[i] - fpr.zl.60.test.10perc[i-1] ) * 
      ( tpr.zl.60.test.10perc[i-1] + (tpr.zl.60.test.10perc[i] - tpr.zl.60.test.10perc[i-1])/2 )
    AUC.temp = sum(AUC.temp, S)
  }
  AUC.zl.60.test.10perc <- c(AUC.zl.60.test.10perc, AUC.temp)
  
  ####### 10.4.4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.60.test.10perc[fpr.zl.60.test.10perc<=0.01])
  x2 <- min(fpr.zl.60.test.10perc[fpr.zl.60.test.10perc> 0.01])
  
  y1 <- tpr.zl.60.test.10perc[fpr.zl.60.test.10perc==x1]
  y2 <- tpr.zl.60.test.10perc[fpr.zl.60.test.10perc==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.10perc <- c(tpr.hat.10perc, tpr.hat.temp)
  
}

se.AUC.zl.60.test.10perc <- sd(AUC.zl.60.test.10perc)/10
AUC.zl.60.test.10perc <- mean(AUC.zl.60.test.10perc)

se.tpr.hat.10perc <- sd(tpr.hat.10perc)/10
tpr.hat.10perc <- mean(tpr.hat.10perc)

################################### 10.6 15%
AUC.zl.60.test.15perc <- c()
tpr.hat.15perc <- c()

for(j in 1:100){
  ####### 10.4.1 组合样本
  set.seed(j)
  temp.60o.p <- test.result.zl.p.60o[sample(1:length(test.result.zl.p.60o),74,replace=FALSE)]
  test.result.zl.60.p.15perc <- c(test.result.zl.p.60y, temp.60o.p)
  
  set.seed(j)
  temp.60y.h <- test.result.zl.h.60y[sample(1:length(test.result.zl.h.60y),11781,replace=FALSE)]
  test.result.zl.60.h.15perc <- c(test.result.zl.h.60o, temp.60y.h)
  
  ####### 10.4.2 fpr,tpr
  tpr.zl.60.test.15perc <- c(0)
  fpr.zl.60.test.15perc <- c(0)
  
  test.result.zl.60.15perc = unique(c(test.result.zl.60.p.15perc, test.result.zl.60.h.15perc))
  
  for(i in sort(test.result.zl.60.15perc, decreasing = T)){
    tpr.zl.60.test.15perc = c(tpr.zl.60.test.15perc, 
                              sum(test.result.zl.60.p.15perc >= i)/length(test.result.zl.60.p.15perc))
    fpr.zl.60.test.15perc = c(fpr.zl.60.test.15perc,
                              sum(test.result.zl.60.h.15perc >= i)/length(test.result.zl.60.h.15perc))
  }
  
  ####### 10.4.3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.60.test.15perc)){
    S = ( fpr.zl.60.test.15perc[i] - fpr.zl.60.test.15perc[i-1] ) * 
      ( tpr.zl.60.test.15perc[i-1] + (tpr.zl.60.test.15perc[i] - tpr.zl.60.test.15perc[i-1])/2 )
    AUC.temp = sum(AUC.temp, S)
  }
  AUC.zl.60.test.15perc <- c(AUC.zl.60.test.15perc, AUC.temp)
  
  ####### 10.4.4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.60.test.15perc[fpr.zl.60.test.15perc<=0.01])
  x2 <- min(fpr.zl.60.test.15perc[fpr.zl.60.test.15perc> 0.01])
  
  y1 <- tpr.zl.60.test.15perc[fpr.zl.60.test.15perc==x1]
  y2 <- tpr.zl.60.test.15perc[fpr.zl.60.test.15perc==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.15perc <- c(tpr.hat.15perc, tpr.hat.temp)
  
}

se.AUC.zl.60.test.15perc <- sd(AUC.zl.60.test.15perc)/10
AUC.zl.60.test.15perc <- mean(AUC.zl.60.test.15perc)

se.tpr.hat.15perc <- sd(tpr.hat.15perc)/10
tpr.hat.15perc <- mean(tpr.hat.15perc)

################################### 10.7 20%
AUC.zl.60.test.20perc <- c()
tpr.hat.20perc <- c()

for(j in 1:100){
  ####### 10.4.1 组合样本
  set.seed(j)
  temp.60o.p <- test.result.zl.p.60o[sample(1:length(test.result.zl.p.60o),106,replace=FALSE)]
  test.result.zl.60.p.20perc <- c(test.result.zl.p.60y, temp.60o.p)
  
  set.seed(j)
  temp.60y.h <- test.result.zl.h.60y[sample(1:length(test.result.zl.h.60y),8316,replace=FALSE)]
  test.result.zl.60.h.20perc <- c(test.result.zl.h.60o, temp.60y.h)
  
  ####### 10.4.2 fpr,tpr
  tpr.zl.60.test.20perc <- c(0)
  fpr.zl.60.test.20perc <- c(0)
  
  test.result.zl.60.20perc = unique(c(test.result.zl.60.p.20perc, test.result.zl.60.h.20perc))
  
  for(i in sort(test.result.zl.60.20perc, decreasing = T)){
    tpr.zl.60.test.20perc = c(tpr.zl.60.test.20perc, 
                              sum(test.result.zl.60.p.20perc >= i)/length(test.result.zl.60.p.20perc))
    fpr.zl.60.test.20perc = c(fpr.zl.60.test.20perc,
                              sum(test.result.zl.60.h.20perc >= i)/length(test.result.zl.60.h.20perc))
  }
  
  ####### 10.4.3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.60.test.20perc)){
    S = ( fpr.zl.60.test.20perc[i] - fpr.zl.60.test.20perc[i-1] ) * 
      ( tpr.zl.60.test.20perc[i-1] + (tpr.zl.60.test.20perc[i] - tpr.zl.60.test.20perc[i-1])/2 )
    AUC.temp = sum(AUC.temp, S)
  }
  AUC.zl.60.test.20perc <- c(AUC.zl.60.test.20perc, AUC.temp)
  
  ####### 10.4.4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.60.test.20perc[fpr.zl.60.test.20perc<=0.01])
  x2 <- min(fpr.zl.60.test.20perc[fpr.zl.60.test.20perc> 0.01])
  
  y1 <- tpr.zl.60.test.20perc[fpr.zl.60.test.20perc==x1]
  y2 <- tpr.zl.60.test.20perc[fpr.zl.60.test.20perc==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.20perc <- c(tpr.hat.20perc, tpr.hat.temp)
  
}

se.AUC.zl.60.test.20perc <- sd(AUC.zl.60.test.20perc)/10
AUC.zl.60.test.20perc <- mean(AUC.zl.60.test.20perc)

se.tpr.hat.20perc <- sd(tpr.hat.20perc)/10
tpr.hat.20perc <- mean(tpr.hat.20perc)

################################### 10.8 25%
AUC.zl.60.test.25perc <- c()
tpr.hat.25perc <- c()

for(j in 1:100){
  ####### 10.4.1 组合样本
  set.seed(j)
  temp.60o.p <- test.result.zl.p.60o[sample(1:length(test.result.zl.p.60o),141,replace=FALSE)]
  test.result.zl.60.p.25perc <- c(test.result.zl.p.60y, temp.60o.p)
  
  set.seed(j)
  temp.60y.h <- test.result.zl.h.60y[sample(1:length(test.result.zl.h.60y),6237,replace=FALSE)]
  test.result.zl.60.h.25perc <- c(test.result.zl.h.60o, temp.60y.h)
  
  ####### 10.4.2 fpr,tpr
  tpr.zl.60.test.25perc <- c(0)
  fpr.zl.60.test.25perc <- c(0)
  
  test.result.zl.60.25perc = unique(c(test.result.zl.60.p.25perc, test.result.zl.60.h.25perc))
  
  for(i in sort(test.result.zl.60.25perc, decreasing = T)){
    tpr.zl.60.test.25perc = c(tpr.zl.60.test.25perc, 
                              sum(test.result.zl.60.p.25perc >= i)/length(test.result.zl.60.p.25perc))
    fpr.zl.60.test.25perc = c(fpr.zl.60.test.25perc,
                              sum(test.result.zl.60.h.25perc >= i)/length(test.result.zl.60.h.25perc))
  }
  
  ####### 10.4.3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.60.test.25perc)){
    S = ( fpr.zl.60.test.25perc[i] - fpr.zl.60.test.25perc[i-1] ) * 
      ( tpr.zl.60.test.25perc[i-1] + (tpr.zl.60.test.25perc[i] - tpr.zl.60.test.25perc[i-1])/2 )
    AUC.temp = sum(AUC.temp, S)
  }
  AUC.zl.60.test.25perc <- c(AUC.zl.60.test.25perc, AUC.temp)
  
  ####### 10.4.4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.60.test.25perc[fpr.zl.60.test.25perc<=0.01])
  x2 <- min(fpr.zl.60.test.25perc[fpr.zl.60.test.25perc> 0.01])
  
  y1 <- tpr.zl.60.test.25perc[fpr.zl.60.test.25perc==x1]
  y2 <- tpr.zl.60.test.25perc[fpr.zl.60.test.25perc==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.25perc <- c(tpr.hat.25perc, tpr.hat.temp)
  
}

se.AUC.zl.60.test.25perc <- sd(AUC.zl.60.test.25perc)/10
AUC.zl.60.test.25perc <- mean(AUC.zl.60.test.25perc)

se.tpr.hat.25perc <- sd(tpr.hat.25perc)/10
tpr.hat.25perc <- mean(tpr.hat.25perc)


################################### 10.9 30%
AUC.zl.60.test.30perc <- c()
tpr.hat.30perc <- c()

for(j in 1:100){
  ####### 10.4.1 组合样本
  set.seed(j)
  temp.60o.p <- test.result.zl.p.60o[sample(1:length(test.result.zl.p.60o),181,replace=FALSE)]
  test.result.zl.60.p.30perc <- c(test.result.zl.p.60y, temp.60o.p)
  
  set.seed(j)
  temp.60y.h <- test.result.zl.h.60y[sample(1:length(test.result.zl.h.60y),4851,replace=FALSE)]
  test.result.zl.60.h.30perc <- c(test.result.zl.h.60o, temp.60y.h)
  
  ####### 10.4.2 fpr,tpr
  tpr.zl.60.test.30perc <- c(0)
  fpr.zl.60.test.30perc <- c(0)
  
  test.result.zl.60.30perc = unique(c(test.result.zl.60.p.30perc, test.result.zl.60.h.30perc))
  
  for(i in sort(test.result.zl.60.30perc, decreasing = T)){
    tpr.zl.60.test.30perc = c(tpr.zl.60.test.30perc, 
                              sum(test.result.zl.60.p.30perc >= i)/length(test.result.zl.60.p.30perc))
    fpr.zl.60.test.30perc = c(fpr.zl.60.test.30perc,
                              sum(test.result.zl.60.h.30perc >= i)/length(test.result.zl.60.h.30perc))
  }
  
  ####### 10.4.3 计算AUC
  AUC.temp <- c()
  for(i in 2:length(fpr.zl.60.test.30perc)){
    S = ( fpr.zl.60.test.30perc[i] - fpr.zl.60.test.30perc[i-1] ) * 
      ( tpr.zl.60.test.30perc[i-1] + (tpr.zl.60.test.30perc[i] - tpr.zl.60.test.30perc[i-1])/2 )
    AUC.temp = sum(AUC.temp, S)
  }
  AUC.zl.60.test.30perc <- c(AUC.zl.60.test.30perc, AUC.temp)
  
  ####### 10.4.4 计算tpr(fpr=1%)
  x1 <- max(fpr.zl.60.test.30perc[fpr.zl.60.test.30perc<=0.01])
  x2 <- min(fpr.zl.60.test.30perc[fpr.zl.60.test.30perc> 0.01])
  
  y1 <- tpr.zl.60.test.30perc[fpr.zl.60.test.30perc==x1]
  y2 <- tpr.zl.60.test.30perc[fpr.zl.60.test.30perc==x2]
  
  tpr.hat.temp <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)
  tpr.hat.30perc <- c(tpr.hat.30perc, tpr.hat.temp)
  
}

se.AUC.zl.60.test.30perc <- sd(AUC.zl.60.test.30perc)/10
AUC.zl.60.test.30perc <- mean(AUC.zl.60.test.30perc)

se.tpr.hat.30perc <- sd(tpr.hat.30perc)/10
tpr.hat.30perc <- mean(tpr.hat.30perc)

################################### 10.10 整理AUC和tpr
AUC.zl.60.test.perc <- c(AUC.zl.60.test.0perc, AUC.zl.60.test.5perc, AUC.zl.60.test.10perc,
                         AUC.zl.60.test.15perc, AUC.zl.60.test.20perc, AUC.zl.60.test.25perc,
                         AUC.zl.60.test.30perc)

se.AUC.zl.60.test.perc <- c(se.AUC.zl.60.test.0perc, se.AUC.zl.60.test.5perc, se.AUC.zl.60.test.10perc,
                         se.AUC.zl.60.test.15perc, se.AUC.zl.60.test.20perc, se.AUC.zl.60.test.25perc,
                         se.AUC.zl.60.test.30perc)

errorbar.perc.AUC.low <- AUC.zl.60.test.perc - se.AUC.zl.60.test.perc
errorbar.perc.AUC.high <- AUC.zl.60.test.perc + se.AUC.zl.60.test.perc


tpr.hat.perc <- c(tpr.hat.0perc, tpr.hat.5perc, tpr.hat.10perc, 
                  tpr.hat.15perc, tpr.hat.20perc, tpr.hat.25perc, 
                  tpr.hat.30perc)

se.tpr.hat.perc <- c(se.tpr.hat.0perc, se.tpr.hat.5perc, se.tpr.hat.10perc, 
                     se.tpr.hat.15perc, se.tpr.hat.20perc, se.tpr.hat.25perc, 
                     se.tpr.hat.30perc)

errorbar.perc.tpr.low <- tpr.hat.perc - se.tpr.hat.perc
errorbar.perc.tpr.high <- tpr.hat.perc + se.tpr.hat.perc


#####

#############################################################
##### (12) CART test (gender)
##########################################################
################################ male
################################ 11.1.1  取出male
data.zl.p.test.male <- data.zl.p.test[data.zl.p.test$sex =="男", c("age", "Alb", "HCT", "LYMPH%", "class")]
data.zl.h.test.male <- data.zl.h.test[data.zl.h.test$sex =="男", c("age", "Alb", "HCT", "LYMPH%", "class")]

################################ 11.1.2 make predictions
test.result.zl.p.male <- predict(fit.zl.CART, data.zl.p.test.male[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.zl.h.male <- predict(fit.zl.CART, data.zl.h.test.male[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]


################################### 11.1.3 计算tpr,fpr

tpr.zl.male.test = c(0)
fpr.zl.male.test = c(0)

test.result.zl.male = unique(c(test.result.zl.p.male, test.result.zl.h.male))

for(i in sort(test.result.zl.male, decreasing = T)){
  tpr.zl.male.test = c(tpr.zl.male.test, sum(test.result.zl.p.male >= i)/length(test.result.zl.p.male))
  fpr.zl.male.test = c(fpr.zl.male.test, sum(test.result.zl.h.male >= i)/length(test.result.zl.h.male))
}

################################## 11.1.4 计算AUC
AUC.zl.male.test = c()
for(i in 2:length(fpr.zl.male.test)){
  S = (fpr.zl.male.test[i]-fpr.zl.male.test[i-1]) * (tpr.zl.male.test[i-1] + (tpr.zl.male.test[i]-tpr.zl.male.test[i-1])/2 )
  AUC.zl.male.test = sum(AUC.zl.male.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.male.test[fpr.zl.male.test<=0.01] )
x2 <- min( fpr.zl.male.test[fpr.zl.male.test> 0.01] )

y1 <- max( tpr.zl.male.test[fpr.zl.male.test==x1] )
y2 <- min( tpr.zl.male.test[fpr.zl.male.test==x2] )

tpr.hat.male.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

################################ female
################################ 11.2.1  取出female
data.zl.p.test.female   = data.zl.p.test[data.zl.p.test$sex =="女", c("age", "Alb", "HCT", "LYMPH%", "class")]
data.zl.h.test.female   = data.zl.h.test[data.zl.h.test$sex =="女", c("age", "Alb", "HCT", "LYMPH%", "class")]

################################ 11.2.2 make predictions
test.result.zl.p.female <- predict(fit.zl.CART, data.zl.p.test.female[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.zl.h.female <- predict(fit.zl.CART, data.zl.h.test.female[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]


################################### 11.2.3 计算tpr,fpr

tpr.zl.female.test = c(0)
fpr.zl.female.test = c(0)

test.result.zl.female = unique(c(test.result.zl.p.female, test.result.zl.h.female))

for(i in sort(test.result.zl.female, decreasing = T)){
  tpr.zl.female.test = c(tpr.zl.female.test, sum(test.result.zl.p.female >= i)/length(test.result.zl.p.female))
  fpr.zl.female.test = c(fpr.zl.female.test, sum(test.result.zl.h.female >= i)/length(test.result.zl.h.female))
}

################################## 11.2.4 计算AUC
AUC.zl.female.test = c()
for(i in 2:length(fpr.zl.female.test)){
  S = (fpr.zl.female.test[i]-fpr.zl.female.test[i-1]) * (tpr.zl.female.test[i-1] + (tpr.zl.female.test[i]-tpr.zl.female.test[i-1])/2 )
  AUC.zl.female.test = sum(AUC.zl.female.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.female.test[fpr.zl.female.test<=0.01] )
x2 <- min( fpr.zl.female.test[fpr.zl.female.test> 0.01] )

y1 <- max( tpr.zl.female.test[fpr.zl.female.test==x1] )
y2 <- min( tpr.zl.female.test[fpr.zl.female.test==x2] )

tpr.hat.female.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)


#####

#############################################################
##### (13) CART test (stage)
##########################################################
################################ early
################################ 11.1.1  取出earlyⅡⅠ期"
data.zl.p.test.early1   = data.zl.p.test[data.zl.p.test$stage =="Ⅰ期", c( "age", "Alb", "HCT", "LYMPH%", "class")]
data.zl.p.test.early2   = data.zl.p.test[data.zl.p.test$stage =="Ⅱ期", c("age", "Alb", "HCT", "LYMPH%", "class")]
data.zl.p.test.early <- rbind(data.zl.p.test.early1, data.zl.p.test.early2)

################################ 11.1.2 make predictions
test.result.zl.p.early <- predict(fit.zl.CART, data.zl.p.test.early[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
### test.result.zl.h.CART

################################### 11.1.3 计算tpr,fpr

tpr.zl.early.test = c(0)
fpr.zl.early.test = c(0)

test.result.zl.early = unique(c(test.result.zl.p.early, test.result.zl.h.CART))

for(i in sort(test.result.zl.early, decreasing = T)){
  tpr.zl.early.test = c(tpr.zl.early.test, sum(test.result.zl.p.early >= i)/length(test.result.zl.p.early))
  fpr.zl.early.test = c(fpr.zl.early.test, sum(test.result.zl.h.CART >= i)/length(test.result.zl.h.CART))
}

################################## 11.1.4 计算AUC
AUC.zl.early.test = c()
for(i in 2:length(fpr.zl.early.test)){
  S = (fpr.zl.early.test[i]-fpr.zl.early.test[i-1]) * (tpr.zl.early.test[i-1] + (tpr.zl.early.test[i]-tpr.zl.early.test[i-1])/2 )
  AUC.zl.early.test = sum(AUC.zl.early.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.early.test[fpr.zl.early.test<=0.01] )
x2 <- min( fpr.zl.early.test[fpr.zl.early.test> 0.01] )

y1 <- max( tpr.zl.early.test[fpr.zl.early.test==x1] )
y2 <- min( tpr.zl.early.test[fpr.zl.early.test==x2] )

tpr.hat.early.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)


################################ advanced
################################ 11.2.1  取出 advanced Ⅲ期,Ⅳ期
data.zl.p.test.advanced3   = data.zl.p.test[data.zl.p.test$stage =="Ⅲ期", c( "age", "Alb", "HCT", "LYMPH%", "class")]
data.zl.p.test.advanced4   = data.zl.p.test[data.zl.p.test$stage =="Ⅳ期", c("age", "Alb", "HCT", "LYMPH%", "class")]
data.zl.p.test.advanced <- rbind(data.zl.p.test.advanced3, data.zl.p.test.advanced4)

################################ 11.2.2 make predictions
test.result.zl.p.advanced <- predict(fit.zl.CART, data.zl.p.test.advanced[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
### test.result.zl.h.CART

################################### 11.2.3 计算tpr,fpr

tpr.zl.advanced.test = c(0)
fpr.zl.advanced.test = c(0)

test.result.zl.advanced = unique(c(test.result.zl.p.advanced, test.result.zl.h.CART))

for(i in sort(test.result.zl.advanced, decreasing = T)){
  tpr.zl.advanced.test = c(tpr.zl.advanced.test, sum(test.result.zl.p.advanced >= i)/length(test.result.zl.p.advanced))
  fpr.zl.advanced.test = c(fpr.zl.advanced.test, sum(test.result.zl.h.CART >= i)/length(test.result.zl.h.CART))
}

################################## 11.2.4 计算AUC
AUC.zl.advanced.test = c()
for(i in 2:length(fpr.zl.advanced.test)){
  S = (fpr.zl.advanced.test[i]-fpr.zl.advanced.test[i-1]) * (tpr.zl.advanced.test[i-1] + (tpr.zl.advanced.test[i]-tpr.zl.advanced.test[i-1])/2 )
  AUC.zl.advanced.test = sum(AUC.zl.advanced.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.advanced.test[fpr.zl.advanced.test<=0.01] )
x2 <- min( fpr.zl.advanced.test[fpr.zl.advanced.test> 0.01] )

y1 <- max( tpr.zl.advanced.test[fpr.zl.advanced.test==x1] )
y2 <- min( tpr.zl.advanced.test[fpr.zl.advanced.test==x2] )

tpr.hat.advanced.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)


################################ Unspecified stages
################################ 11.3.1  取出Unspecified stages
data.zl.p.test.Unspecified   = data.zl.p.test[data.zl.p.test$stage =="unknown", c( "age", "Alb", "HCT", "LYMPH%", "class")]

################################ 11.3.2 make predictions
test.result.zl.p.Unspecified <- predict(fit.zl.CART, data.zl.p.test.Unspecified[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
### test.result.zl.h.CART

################################### 11.3.3 计算tpr,fpr

tpr.zl.Unspecified.test = c(0)
fpr.zl.Unspecified.test = c(0)

test.result.zl.Unspecified = unique(c(test.result.zl.p.Unspecified, test.result.zl.h.CART))

for(i in sort(test.result.zl.Unspecified, decreasing = T)){
  tpr.zl.Unspecified.test = c(tpr.zl.Unspecified.test, sum(test.result.zl.p.Unspecified >= i)/length(test.result.zl.p.Unspecified))
  fpr.zl.Unspecified.test = c(fpr.zl.Unspecified.test, sum(test.result.zl.h.CART >= i)/length(test.result.zl.h.CART))
}

################################## 11.3.4 计算AUC
AUC.zl.Unspecified.test = c()
for(i in 2:length(fpr.zl.Unspecified.test)){
  S = (fpr.zl.Unspecified.test[i]-fpr.zl.Unspecified.test[i-1]) * (tpr.zl.Unspecified.test[i-1] + (tpr.zl.Unspecified.test[i]-tpr.zl.Unspecified.test[i-1])/2 )
  AUC.zl.Unspecified.test = sum(AUC.zl.Unspecified.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.Unspecified.test[fpr.zl.Unspecified.test<=0.01] )
x2 <- min( fpr.zl.Unspecified.test[fpr.zl.Unspecified.test> 0.01] )

y1 <- max( tpr.zl.Unspecified.test[fpr.zl.Unspecified.test==x1] )
y2 <- min( tpr.zl.Unspecified.test[fpr.zl.Unspecified.test==x2] )

tpr.hat.Unspecified.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (14) CART test (location)
##########################################################
################################ colon
################################ 14.1.1  取出colon
data.zl.p.test.colon   = data.zl.p.test[data.zl.p.test$location =="colon", c("age", "Alb", "HCT", "LYMPH%", "class")]

################################ 14.1.2 make predictions
test.result.zl.p.colon <- predict(fit.zl.CART, data.zl.p.test.colon[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]


################################### 14.1.3 计算tpr,fpr

tpr.zl.colon.test = c(0)
fpr.zl.colon.test = c(0)

test.result.zl.colon = unique(c(test.result.zl.p.colon, test.result.zl.h.CART))

for(i in sort(test.result.zl.colon, decreasing = T)){
  tpr.zl.colon.test = c(tpr.zl.colon.test, sum(test.result.zl.p.colon >= i)/length(test.result.zl.p.colon))
  fpr.zl.colon.test = c(fpr.zl.colon.test, sum(test.result.zl.h.CART >= i)/length(test.result.zl.h.CART))
}

################################## 14.1.4 计算AUC
AUC.zl.colon.test = c()
for(i in 2:length(fpr.zl.colon.test)){
  S = (fpr.zl.colon.test[i]-fpr.zl.colon.test[i-1]) * (tpr.zl.colon.test[i-1] + (tpr.zl.colon.test[i]-tpr.zl.colon.test[i-1])/2 )
  AUC.zl.colon.test = sum(AUC.zl.colon.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.colon.test[fpr.zl.colon.test<=0.01] )
x2 <- min( fpr.zl.colon.test[fpr.zl.colon.test> 0.01] )

y1 <- max( tpr.zl.colon.test[fpr.zl.colon.test==x1] )
y2 <- min( tpr.zl.colon.test[fpr.zl.colon.test==x2] )

tpr.hat.colon.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

################################ rectum
################################ 14.2.1  取出rectum
data.zl.p.test.rectum   = data.zl.p.test[data.zl.p.test$location =="rectum", c("age", "Alb", "HCT", "LYMPH%", "class")]

################################ 14.2.2 make predictions
test.result.zl.p.rectum <- predict(fit.zl.CART, data.zl.p.test.rectum[, c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]


################################### 14.2.3 计算tpr,fpr

tpr.zl.rectum.test = c(0)
fpr.zl.rectum.test = c(0)

test.result.zl.rectum = unique(c(test.result.zl.p.rectum, test.result.zl.h.CART))

for(i in sort(test.result.zl.rectum, decreasing = T)){
  tpr.zl.rectum.test = c(tpr.zl.rectum.test, sum(test.result.zl.p.rectum >= i)/length(test.result.zl.p.rectum))
  fpr.zl.rectum.test = c(fpr.zl.rectum.test, sum(test.result.zl.h.CART >= i)/length(test.result.zl.h.CART))
}

################################## 14.2.4 计算AUC
AUC.zl.rectum.test = c()
for(i in 2:length(fpr.zl.rectum.test)){
  S = (fpr.zl.rectum.test[i]-fpr.zl.rectum.test[i-1]) * (tpr.zl.rectum.test[i-1] + (tpr.zl.rectum.test[i]-tpr.zl.rectum.test[i-1])/2 )
  AUC.zl.rectum.test = sum(AUC.zl.rectum.test, S)
}

### tpr(fpr=1%)
x1 <- max( fpr.zl.rectum.test[fpr.zl.rectum.test<=0.01] )
x2 <- min( fpr.zl.rectum.test[fpr.zl.rectum.test> 0.01] )

y1 <- max( tpr.zl.rectum.test[fpr.zl.rectum.test==x1] )
y2 <- min( tpr.zl.rectum.test[fpr.zl.rectum.test==x2] )

tpr.hat.rectum.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (15) sg
##########################################################
#####################  预处理
head(patient.info.sg)
colnames(patient.info.sg)
patient.matrix.sg = cbind( patient.info.sg[, c(5,24,27,20)], 
                           class = rep( "patient", dim(patient.info.sg)[1] ),
                           stringsAsFactors = FALSE
)
colnames(patient.matrix.sg) = c("age", "Alb", "HCT", "LYMPH%", "class")
head(patient.matrix.sg)

head(health.info.sg)
colnames(health.info.sg)
health.matrix.sg  = cbind( health.info.sg[, c(3,22,25,18)], 
                           class = rep( "health", dim(health.info.sg)[1] ),
                           stringsAsFactors = FALSE
                           )
colnames(health.matrix.sg) = c("age", "Alb", "HCT", "LYMPH%", "class")
head(health.matrix.sg)

colnames(patient.matrix.sg)
colnames(health.matrix.sg)

#################### 数据分割
dim.sg.p = dim(patient.matrix.sg)[1]
dim.sg.h = dim(health.matrix.sg)[1]

dim.sg.p.test = floor(dim.sg.p*0.3)
dim.sg.h.test = floor(dim.sg.h*0.3)

set.seed(100)
data.label.sg.p.test  = sample( seq(1,dim.sg.p), dim.sg.p.test, replace=FALSE )
data.label.sg.p.train = setdiff( seq(1,dim.sg.p), data.label.sg.p.test )

set.seed(100)
data.label.sg.h.test  = sample( seq(1,dim.sg.h), dim.sg.h.test, replace=FALSE )
data.label.sg.h.train = setdiff( seq(1,dim.sg.h), data.label.sg.h.test )

data.sg.p.test = patient.matrix.sg[data.label.sg.p.test, ]
data.sg.h.test = health.matrix.sg[data.label.sg.h.test, ]

data.sg.p.train = patient.matrix.sg[data.label.sg.p.train, ]
data.sg.h.train = health.matrix.sg[data.label.sg.h.train, ]
data.sg.train   = rbind( data.sg.p.train, data.sg.h.train )

####################### fit.zl.CART, sgtest
test.result.zl.sg.p.CART <- predict(fit.zl.CART, data.sg.p.test[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.zl.sg.h.CART <- predict(fit.zl.CART, data.sg.h.test[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]

tpr.zl.sg.test = c(0)
fpr.zl.sg.test = c(0)

test.result.zl.sg.CART = unique(c(test.result.zl.sg.p.CART, test.result.zl.sg.h.CART))
for(i in sort(test.result.zl.sg.CART,decreasing=TRUE)){
  tpr.zl.sg.test = c(tpr.zl.sg.test, sum(test.result.zl.sg.p.CART >= i)/length(test.result.zl.sg.p.CART))
  fpr.zl.sg.test = c(fpr.zl.sg.test, sum(test.result.zl.sg.h.CART >= i)/length(test.result.zl.sg.h.CART))
}

AUC.zl.sg.CART.test = c()
for(i in 2:length(fpr.zl.sg.test)){
  S <- ( fpr.zl.sg.test[i]-fpr.zl.sg.test[i-1] ) * 
    (tpr.zl.sg.test[i-1] + (tpr.zl.sg.test[i]-tpr.zl.sg.test[i-1])/2 )
  AUC.zl.sg.CART.test = sum(AUC.zl.sg.CART.test, S)
}


########################## fit.CART.sg, sgtest
fit.sg.CART = rpart(class~., data=data.sg.train)

test.result.sg.sg.p.CART <- predict(fit.sg.CART, data.sg.p.test[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.sg.sg.h.CART <- predict(fit.sg.CART, data.sg.h.test[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]

tpr.sg.sg.test = c(0)
fpr.sg.sg.test = c(0)

test.result.sg.sg.CART = unique(c(test.result.sg.sg.p.CART, test.result.sg.sg.h.CART))
for(i in sort(test.result.sg.sg.CART,decreasing=TRUE)){
  tpr.sg.sg.test = c(tpr.sg.sg.test, sum(test.result.sg.sg.p.CART >= i)/length(test.result.sg.sg.p.CART))
  fpr.sg.sg.test = c(fpr.sg.sg.test, sum(test.result.sg.sg.h.CART >= i)/length(test.result.sg.sg.h.CART))
}

AUC.sg.sg.CART.test = c()
for(i in 2:length(fpr.sg.sg.test)){
  S <- ( fpr.sg.sg.test[i]-fpr.sg.sg.test[i-1] ) * 
    (tpr.sg.sg.test[i-1] + (tpr.sg.sg.test[i]-tpr.sg.sg.test[i-1])/2 )
  AUC.sg.sg.CART.test = sum(AUC.sg.sg.CART.test, S)
}

x1 <- max( fpr.sg.sg.test[fpr.sg.sg.test<=0.01] )
x2 <- min( fpr.sg.sg.test[fpr.sg.sg.test> 0.01] )

y1 <- max( tpr.sg.sg.test[fpr.sg.sg.test==x1] )
y2 <- min( tpr.sg.sg.test[fpr.sg.sg.test==x2] )

tpr.hat.sg.sg.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)

#####

#############################################################
##### (16) 采样一个更general的population来看age model的AUC
##########################################################
####################################################### 采样
set.seed(100)
data.zl.p.test.general <- rbind(data.zl.p.test.60y, 
                                data.zl.p.test.60o[sample(1:(dim(data.zl.p.test.60o)[1]),150),])
data.zl.h.test.general <- rbind(data.zl.h.test.60y, data.zl.h.test.60o)

###################################################### age model
tpr.zl.age.general <- c( sum(data.zl.p.test.general$age >= 100)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 95)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 90)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 85)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 80)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 75)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 70)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 65)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 60)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 55)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 50)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 45)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 42)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 40)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 35)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 30)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 25)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 20)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 15)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 10)/length(data.zl.p.test.general$age),
                       sum(data.zl.p.test.general$age >= 5 )/length(data.zl.p.test.general$age)
)

fpr.zl.age.general <- c( sum(data.zl.h.test.general$age >= 100)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 95)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 90)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 85)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 80)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 75)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 70)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 65)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 60)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 55)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 50)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 45)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 42)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 40)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 35)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 30)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 25)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 20)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 15)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 10)/length(data.zl.h.test.general$age),
                       sum(data.zl.h.test.general$age >= 5 )/length(data.zl.h.test.general$age)
)

#### AUC
AUC.zl.age.general <- c()
for(i in 2:length(fpr.zl.age.general)){
  S <-( fpr.zl.age.general[i]-fpr.zl.age.general[i-1] ) * 
    ( tpr.zl.age.general[i-1] + (tpr.zl.age.general[i]-tpr.zl.age.general[i-1])/2 )
  AUC.zl.age.general <- sum(AUC.zl.age.general, S)
}

################################################################ CART model
test.result.zl.p.CART.general <- predict(fit.zl.CART, data.zl.p.test.general, type="prob")[,2]
test.result.zl.h.CART.general <- test.result.zl.h.CART

################################### 9.2 计算tpr,fpr
tpr.zl.CART.general <- c(0)
fpr.zl.CART.general <- c(0)

test.result.zl.CART.general <- unique(c(test.result.zl.p.CART.general, test.result.zl.h.CART.general))

for(i in sort(test.result.zl.CART.general, decreasing=TRUE)){
  tpr.zl.CART.general <- c(tpr.zl.CART.general, sum(test.result.zl.p.CART.general >= i)/length(test.result.zl.p.CART.general))
  fpr.zl.CART.general <- c(fpr.zl.CART.general, sum(test.result.zl.h.CART.general >= i)/length(test.result.zl.h.CART.general))
}

################################## 9.3 计算AUC
AUC.zl.CART.general <- c()
for(i in 2:length(fpr.zl.CART.general)){
  S <- ( fpr.zl.CART.general[i]-fpr.zl.CART.general[i-1] ) * 
    ( tpr.zl.CART.general[i-1] + (tpr.zl.CART.general[i]-tpr.zl.CART.general[i-1])/2 )
  AUC.zl.CART.general <- sum(AUC.zl.CART.general, S)
}


#####

#############################################################
##### (17) 用结肠1-3期和直肠1期训练和测试 
##########################################################
###################### 16.1 提取训练集
index1 <- (data.zl.p.train$location=="colon") &
  ((data.zl.p.train$stage=="Ⅰ期")|(data.zl.p.train$stage=="Ⅱ期")|(data.zl.p.train$stage=="Ⅲ期"))
# index1 <- (data.zl.p.train$location=="colon") &
#   ((data.zl.p.train$stage=="Ⅰ期")|(data.zl.p.train$stage=="Ⅱ期"))
index2 <- (data.zl.p.train$location=="rectum") & (data.zl.p.train$stage=="Ⅰ期")
index12 <- (index1 | index2)
# index12 <- index1
data.zl.p.train.colon123.rec1 <- data.zl.p.train[index12, c("age", "Alb", "HCT", "LYMPH%", "class")]##-c(1,2,3)
data.zl.train.colon123.rec1 <- rbind(data.zl.p.train.colon123.rec1, 
                                     data.zl.h.train[,c("age", "Alb", "HCT", "LYMPH%", "class")])
# set.seed(300)
# data.zl.train.colon123.rec1 <- rbind(data.zl.p.train.colon123.rec1, 
#                                      data.zl.h.train[sample(1:(dim(data.zl.h.train)[1]),10000),c("age", "Alb", "HCT", "LYMPH%", "class")])

###################### 16.2 训练
fit.zl.CART.colon123.rec1 <- rpart(class~., data=data.zl.train.colon123.rec1)
summary(fit.zl.CART.colon123.rec1)

#################### 16.3 提取测试集
index3 <- (data.zl.p.test$location=="colon") &
  ((data.zl.p.test$stage=="Ⅰ期")|(data.zl.p.test$stage=="Ⅱ期")|(data.zl.p.test$stage=="Ⅲ期"))
# index3 <- (data.zl.p.test$location=="colon") &
#   ((data.zl.p.test$stage=="Ⅰ期")|(data.zl.p.test$stage=="Ⅱ期"))

index4 <- (data.zl.p.test$location=="rectum") & (data.zl.p.test$stage=="Ⅰ期")
index34 <- (index3 | index4)
# index34 <- index3

data.zl.p.test.colon123.rec1 <- data.zl.p.test[index34, c("age", "Alb", "HCT", "LYMPH%", "class")]##-c(1,2,3)
data.zl.h.test.colon123.rec1 <- data.zl.h.test

#################### 16.3 预测 colon123.rec1
test.result.zl.p.CART.colon123.rec1 <-
  predict(fit.zl.CART.colon123.rec1, data.zl.p.test.colon123.rec1[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.zl.h.CART.colon123.rec1 <-
  predict(fit.zl.CART.colon123.rec1, data.zl.h.test.colon123.rec1[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]

test.result.zl.p.CART.colon123.rec1 <- 
  predict(fit.zl.CART.colon123.rec1, data.zl.p.test[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]
test.result.zl.h.CART.colon123.rec1 <- 
  predict(fit.zl.CART.colon123.rec1, data.zl.h.test[,c("age", "Alb", "HCT", "LYMPH%")], type="prob")[,2]


################################### 16.4 计算tpr,fpr

tpr.zl.CART.colon123.rec1.test = c(0)
fpr.zl.CART.colon123.rec1.test = c(0)

test.result.zl.CART.colon123.rec1 = unique(c(test.result.zl.p.CART.colon123.rec1, test.result.zl.h.CART.colon123.rec1))

for(i in sort(test.result.zl.CART.colon123.rec1, decreasing = T)){
  tpr.zl.CART.colon123.rec1.test = c(tpr.zl.CART.colon123.rec1.test, sum(test.result.zl.p.CART.colon123.rec1 >= i)/length(test.result.zl.p.CART.colon123.rec1))
  fpr.zl.CART.colon123.rec1.test = c(fpr.zl.CART.colon123.rec1.test, sum(test.result.zl.h.CART.colon123.rec1 >= i)/length(test.result.zl.h.CART.colon123.rec1))
}

################################## 16.5 计算AUC
AUC.zl.CART.colon123.rec1.test = c()
for(i in 2:length(fpr.zl.CART.colon123.rec1.test)){
  S = (fpr.zl.CART.colon123.rec1.test[i]-fpr.zl.CART.colon123.rec1.test[i-1]) * (tpr.zl.CART.colon123.rec1.test[i-1] + (tpr.zl.CART.colon123.rec1.test[i]-tpr.zl.CART.colon123.rec1.test[i-1])/2 )
  AUC.zl.CART.colon123.rec1.test = sum(AUC.zl.CART.colon123.rec1.test, S)
}

################################## 16.6 计算tpr(fpr=1%)
x1 <- max( fpr.zl.CART.colon123.rec1.test[fpr.zl.CART.colon123.rec1.test<=0.01] )
x2 <- min( fpr.zl.CART.colon123.rec1.test[fpr.zl.CART.colon123.rec1.test> 0.01] )

y1 <- max( tpr.zl.CART.colon123.rec1.test[fpr.zl.CART.colon123.rec1.test==x1] )
y2 <- min( tpr.zl.CART.colon123.rec1.test[fpr.zl.CART.colon123.rec1.test==x2] )

tpr.hat.colon123.rec1.test <- y1 + (0.01-x1)/(x2-x1)*(y2-y1)


###########################################################
save.image(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/analysis/zl_hospital_analysis_0623.RData"
  )

