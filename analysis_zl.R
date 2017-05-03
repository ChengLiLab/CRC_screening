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
  "/lustre/user/liclab/shiq/CRC/zl_hospital/zl_hospital_preprocess_result.RData"
  )

#####

#############################################################
##### (1) 预处理加类名
##########################################################
head(patient.info.zl)
patient.matrix.zl = cbind( patient.info.zl[, -c(1,2,3)], 
                             class = rep( "patient", dim(patient.info.zl)[1] ),
                             stringsAsFactors = FALSE
                             )
head(patient.matrix.zl)

head(health.info.zl)
health.matrix.zl  = cbind( health.info.zl[, -c(1)], 
                             class = rep( "healthy", dim(health.info.zl)[1] ),
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
dim.p = dim(patient.matrix.zl)[1]
dim.h = dim(health.matrix.zl )[1]

dim.test.p = floor(dim.p*0.2)
dim.test.h = floor(dim.h*0.2)

set.seed(2000)
data.label.test.p  = sample( seq(1,dim.p), dim.test.p, replace=FALSE )
data.label.train.p = setdiff( seq(1,dim.p), data.label.test.p )

set.seed(2000)
data.label.test.h  = sample( seq(1,dim.h), dim.test.h, replace=FALSE )
data.label.train.h = setdiff( seq(1,dim.h), data.label.test.h )

data.test.p = patient.matrix.zl[data.label.test.p, ]
data.test.h = health.matrix.zl[data.label.test.h, ]
data.test   = rbind( data.test.p, data.test.h )

data.train.p = patient.matrix.zl[data.label.train.p, ]
data.train.h = health.matrix.zl[data.label.train.h, ]
data.train   = rbind( data.train.p, data.train.h )

dim.train.p = dim.p-dim.test.p
dim.train.h = dim.h-dim.test.h

dim.valid.p = floor(dim.train.p*0.25)
dim.valid.h = floor(dim.train.h*0.25)


set.seed(3000)
data.label.valid.p  = sample( seq(1,dim.train.p), dim.valid.p, replace=FALSE )
data.label.train.p = setdiff( seq(1,dim.train.p), data.label.valid.p )

set.seed(3000)
data.label.valid.h  = sample( seq(1,dim.train.h), dim.valid.h, replace=FALSE )
data.label.train.h = setdiff( seq(1,dim.train.h), data.label.valid.h )

data.valid.p = data.train.p[data.label.valid.p, ]
data.valid.h = data.train.h[data.label.valid.h, ]
data.valid   = rbind( data.valid.p, data.valid.h )

data.train.p = data.train.p[data.label.train.p, ]
data.train.h = data.train.h[data.label.train.h, ]
data.train   = rbind( data.train.p, data.train.h )



#####

#############################################################
##### (4) age model with training data
##########################################################
tpr.age = c( sum(data.valid.p$age >= 90)/length(data.valid.p$age),
             sum(data.valid.p$age >= 85)/length(data.valid.p$age),
             sum(data.valid.p$age >= 80)/length(data.valid.p$age),
             sum(data.valid.p$age >= 75)/length(data.valid.p$age),
             sum(data.valid.p$age >= 70)/length(data.valid.p$age),
             sum(data.valid.p$age >= 65)/length(data.valid.p$age),
             sum(data.valid.p$age >= 60)/length(data.valid.p$age),
             sum(data.valid.p$age >= 55)/length(data.valid.p$age),
             sum(data.valid.p$age >= 50)/length(data.valid.p$age),
             sum(data.valid.p$age >= 45)/length(data.valid.p$age),
             sum(data.valid.p$age >= 42)/length(data.valid.p$age),
             sum(data.valid.p$age >= 40)/length(data.valid.p$age),
             sum(data.valid.p$age >= 35)/length(data.valid.p$age),
             sum(data.valid.p$age >= 30)/length(data.valid.p$age),
             sum(data.valid.p$age >= 25)/length(data.valid.p$age),
             sum(data.valid.p$age >= 20)/length(data.valid.p$age)
             )

fpr.age = c( sum(data.valid.h$age >= 90)/length(data.valid.h$age),
             sum(data.valid.h$age >= 85)/length(data.valid.h$age),
             sum(data.valid.h$age >= 80)/length(data.valid.h$age),
             sum(data.valid.h$age >= 75)/length(data.valid.h$age),
             sum(data.valid.h$age >= 70)/length(data.valid.h$age),
             sum(data.valid.h$age >= 65)/length(data.valid.h$age),
             sum(data.valid.h$age >= 60)/length(data.valid.h$age),
             sum(data.valid.h$age >= 55)/length(data.valid.h$age),
             sum(data.valid.h$age >= 50)/length(data.valid.h$age),
             sum(data.valid.h$age >= 45)/length(data.valid.h$age),
             sum(data.valid.h$age >= 42)/length(data.valid.h$age),
             sum(data.valid.h$age >= 40)/length(data.valid.h$age),
             sum(data.valid.h$age >= 35)/length(data.valid.h$age),
             sum(data.valid.h$age >= 30)/length(data.valid.h$age),
             sum(data.valid.h$age >= 25)/length(data.valid.h$age),
             sum(data.valid.h$age >= 20)/length(data.valid.h$age)
             )


AUC.age = c()
for(i in 2:length(fpr.age)){
  S   = (fpr.age[i]-fpr.age[i-1]) * (tpr.age[i-1] + (tpr.age[i]-tpr.age[i-1])/2 )
  AUC.age = sum(AUC.age, S)
  }

#####

#############################################################
##### (5) NB training
##########################################################
##################################### 4.1 training with CV
library(e1071)

# fit model
fit.NB <- naiveBayes(class~., data=data.train)

# summarize the fit
summary(fit.NB)

# make predictions
pred.NB <- predict(fit.NB, data.valid[,1:40], type="raw")

table(data.valid$class)

####################################### 4.2 整理预测结果
valid.result.NB.p = pred.NB[1:table(data.valid$class)[2], 2]

valid.result.NB.h = pred.NB[(table(data.valid$class)[2]+1):dim(data.valid)[1], 2]
                  

################################### 4.3 计算tpr,fpr

tpr.NB = c(0)
fpr.NB = c(0)

valid.result.NB = unique(c(valid.result.NB.p, valid.result.NB.h))

for(i in sort(valid.result.NB, decreasing = T)){
  tpr.NB = c(tpr.NB, sum(valid.result.NB.p >= i)/length(valid.result.NB.p))
  fpr.NB = c(fpr.NB, sum(valid.result.NB.h >= i)/length(valid.result.NB.h))
}

################################## 4.4 计算AUC
AUC.NB = c()
for(i in 2:length(fpr.NB)){
  S = (fpr.NB[i]-fpr.NB[i-1]) * (tpr.NB[i-1] + (tpr.NB[i]-tpr.NB[i-1])/2 )
  AUC.NB = sum(AUC.NB, S)
}

#####

#############################################################
##### (6) CART training
##########################################################
library(rpart)

# select best CART

# fit final model
#c(2,4,9,13,18,41)
fit.CART.zl <- rpart(class~., data=data.train[, c("age", "Alb", "CEA", "HCT", "LYMPH%", "class")])
#c(2,4,9,13,18,41)

# summarize the fit
summary(fit.CART.zl)

# make predictions
pred.CART.zl <- predict(fit.CART.zl, data.valid[, c("age", "Alb", "CEA", "HCT", "LYMPH%")], type="prob")

####################################### 4.2 整理预测结果
table(data.valid$class)
valid.result.CART.p.zl = pred.CART.zl[1:table(data.valid$class)[2], 2]
valid.result.CART.h.zl = pred.CART.zl[(table(data.valid$class)[2]+1):dim(data.valid)[1], 2]
                  
################################### 4.3 计算tpr,fpr

tpr.CART.zl = c(0)
fpr.CART.zl = c(0)

valid.result.CART.zl = unique(c(valid.result.CART.p.zl, valid.result.CART.h.zl))

for(i in sort(valid.result.CART.zl, decreasing = T)){
  tpr.CART.zl = c(tpr.CART.zl, sum(valid.result.CART.p.zl >= i)/length(valid.result.CART.p.zl))
  fpr.CART.zl = c(fpr.CART.zl, sum(valid.result.CART.h.zl >= i)/length(valid.result.CART.h.zl))
}

################################## 4.4 计算AUC

AUC.CART.zl = c()
for(i in 2:length(fpr.CART.zl)){
  S = (fpr.CART.zl[i]-fpr.CART.zl[i-1]) * (tpr.CART.zl[i-1] + (tpr.CART.zl[i]-tpr.CART.zl[i-1])/2 )
  AUC.CART.zl = sum(AUC.CART.zl, S)
}


#####

#############################################################
##### (7) age model with testing data
##########################################################
tpr.age.test = c( sum(data.test.p$age >= 90)/length(data.test.p$age),
             sum(data.test.p$age >= 85)/length(data.test.p$age),
             sum(data.test.p$age >= 80)/length(data.test.p$age),
             sum(data.test.p$age >= 75)/length(data.test.p$age),
             sum(data.test.p$age >= 70)/length(data.test.p$age),
             sum(data.test.p$age >= 65)/length(data.test.p$age),
             sum(data.test.p$age >= 60)/length(data.test.p$age),
             sum(data.test.p$age >= 55)/length(data.test.p$age),
             sum(data.test.p$age >= 50)/length(data.test.p$age),
             sum(data.test.p$age >= 45)/length(data.test.p$age),
             sum(data.test.p$age >= 42)/length(data.test.p$age),
             sum(data.test.p$age >= 40)/length(data.test.p$age),
             sum(data.test.p$age >= 35)/length(data.test.p$age),
             sum(data.test.p$age >= 30)/length(data.test.p$age),
             sum(data.test.p$age >= 25)/length(data.test.p$age),
             sum(data.test.p$age >= 20)/length(data.test.p$age)
)

fpr.age.test = c( sum(data.test.h$age >= 90)/length(data.test.h$age),
             sum(data.test.h$age >= 85)/length(data.test.h$age),
             sum(data.test.h$age >= 80)/length(data.test.h$age),
             sum(data.test.h$age >= 75)/length(data.test.h$age),
             sum(data.test.h$age >= 70)/length(data.test.h$age),
             sum(data.test.h$age >= 65)/length(data.test.h$age),
             sum(data.test.h$age >= 60)/length(data.test.h$age),
             sum(data.test.h$age >= 55)/length(data.test.h$age),
             sum(data.test.h$age >= 50)/length(data.test.h$age),
             sum(data.test.h$age >= 45)/length(data.test.h$age),
             sum(data.test.h$age >= 42)/length(data.test.h$age),
             sum(data.test.h$age >= 40)/length(data.test.h$age),
             sum(data.test.h$age >= 35)/length(data.test.h$age),
             sum(data.test.h$age >= 30)/length(data.test.h$age),
             sum(data.test.h$age >= 25)/length(data.test.h$age),
             sum(data.test.h$age >= 20)/length(data.test.h$age)
)


AUC.age.test = c()
for(i in 2:length(fpr.age.test)){
  S = (fpr.age.test[i]-fpr.age.test[i-1]) * (tpr.age.test[i-1] + (tpr.age.test[i]-tpr.age.test[i-1])/2 )
  AUC.age.test = sum(AUC.age.test, S)
}

#####

#############################################################
##### (8) NB testing
##########################################################
##################################### 8.1 predicting

# make predictions
pred.NB.test <- predict(fit.NB, data.test[,1:40], type="raw")

table(data.test$class)

####################################### 8.2 整理预测结果
test.result.NB.p = pred.NB.test[1:table(data.test$class)[2], 2]

test.result.NB.h = pred.NB.test[(table(data.test$class)[2]+1):dim(data.test)[1], 2]


################################### 8.3 计算tpr,fpr

tpr.NB.test = c(0)
fpr.NB.test = c(0)

test.result.NB = unique(c(test.result.NB.p, test.result.NB.h))

for(i in sort(test.result.NB, decreasing = T)){
  tpr.NB.test = c(tpr.NB.test, sum(test.result.NB.p >= i)/length(test.result.NB.p))
  fpr.NB.test = c(fpr.NB.test, sum(test.result.NB.h >= i)/length(test.result.NB.h))
}

################################## 8.4 计算AUC
AUC.NB.test = c()
for(i in 2:length(fpr.NB.test)){
  S = (fpr.NB.test[i]-fpr.NB.test[i-1]) * (tpr.NB.test[i-1] + (tpr.NB.test[i]-tpr.NB.test[i-1])/2 )
  AUC.NB.test = sum(AUC.NB.test, S)
}

#####

#############################################################
##### (9) CART testing
##########################################################

# make predictions
pred.CART.zl.test <- predict(fit.CART.zl, data.test[, c("age", "Alb", "CEA", "HCT", "LYMPH%")], type="prob")

####################################### 4.2 整理预测结果
table(data.test$class)
test.result.CART.p.zl = pred.CART.zl.test[1:table(data.test$class)[2], 2]
test.result.CART.h.zl = pred.CART.zl.test[(table(data.test$class)[2]+1):dim(data.test)[1], 2]

################################### 4.3 计算tpr,fpr

tpr.CART.zl.test = c(0)
fpr.CART.zl.test = c(0)

test.result.CART.zl = unique(c(test.result.CART.p.zl, test.result.CART.h.zl))

for(i in sort(test.result.CART.zl, decreasing = T)){
  tpr.CART.zl.test = c(tpr.CART.zl.test, sum(test.result.CART.p.zl >= i)/length(test.result.CART.p.zl))
  fpr.CART.zl.test = c(fpr.CART.zl.test, sum(test.result.CART.h.zl >= i)/length(test.result.CART.h.zl))
}

################################## 4.4 计算AUC

AUC.CART.zl.test = c()
for(i in 2:length(fpr.CART.zl.test)){
  S = (fpr.CART.zl.test[i]-fpr.CART.zl.test[i-1]) * (tpr.CART.zl.test[i-1] + (tpr.CART.zl.test[i]-tpr.CART.zl.test[i-1])/2 )
  AUC.CART.zl.test = sum(AUC.CART.zl.test, S)
}


#####

#############################################################
##### (10) CART test (age>50 | age<=50)
##########################################################
################################ age>50
data.test.50g   = data.test[data.test$age >50, c("age", "Alb", "CEA", "HCT", "LYMPH%", "class")]
# make predictions
pred.test.50g <- predict(fit.CART.zl, data.test.50g[, -6], type="prob")

########## 整理预测结果
table(data.test.50g$class)
test.result.p.50g = pred.test.50g[1:table(data.test.50g$class)[2], 2]
test.result.h.50g = pred.test.50g[(table(data.test.50g$class)[2]+1):dim(data.test.50g)[1], 2]

########## 计算tpr,fpr
tpr.test.50g = c(0)
fpr.test.50g = c(0)

test.result.50g = unique(c(test.result.p.50g, test.result.h.50g))

for(i in sort(test.result.50g, decreasing = T)){
  tpr.test.50g = c(tpr.test.50g, sum(test.result.p.50g >= i)/length(test.result.p.50g))
  fpr.test.50g = c(fpr.test.50g, sum(test.result.h.50g >= i)/length(test.result.h.50g))
}

########### 计算AUC
AUC.test.50g = c()
for(i in 2:length(fpr.test.50g)){
  S = (fpr.test.50g[i]-fpr.test.50g[i-1]) * (tpr.test.50g[i-1] + (tpr.test.50g[i]-tpr.test.50g[i-1])/2 )
  AUC.test.50g = sum(AUC.test.50g, S)
}


################################ age<=50

data.test.50l   = data.test[data.test$age <= 50, c("age", "Alb", "CEA", "HCT", "LYMPH%", "class")]
# make predictions
pred.test.50l <- predict(fit.CART.zl, data.test.50l[, -6], type="prob")

########## 整理预测结果
table(data.test.50l$class)
test.result.p.50l = pred.test.50l[1:table(data.test.50l$class)[2], 2]
test.result.h.50l = pred.test.50l[(table(data.test.50l$class)[2]+1):dim(data.test.50l)[1], 2]

########## 计算tpr,fpr
tpr.test.50l = c(0)
fpr.test.50l = c(0)

test.result.50l = unique(c(test.result.p.50l, test.result.h.50l))

for(i in sort(test.result.50l, decreasing = T)){
  tpr.test.50l = c(tpr.test.50l, sum(test.result.p.50l >= i)/length(test.result.p.50l))
  fpr.test.50l = c(fpr.test.50l, sum(test.result.h.50l >= i)/length(test.result.h.50l))
}

########### 计算AUC
AUC.test.50l = c()
for(i in 2:length(fpr.test.50l)){
  S = (fpr.test.50l[i]-fpr.test.50l[i-1]) * (tpr.test.50l[i-1] + (tpr.test.50l[i]-tpr.test.50l[i-1])/2 )
  AUC.test.50l = sum(AUC.test.50l, S)
}

#####

#############################################################
##### (10) CART test (gender)
##########################################################
################################ male
data.test.male <- data.test[data.test$sex=="男", c("age", "Alb", "CEA", "HCT", "LYMPH%", "class")]
pred.test.male <- predict(fit.CART.zl, data.test.male[, -6], type="prob")

table(data.test.male$class)
test.result.p.male = pred.test.male[1:table(data.test.male$class)[2],2]
test.result.h.male = pred.test.male[(table(data.test.male$class)[2]+1):dim(data.test.male)[1], 2]

tpr.test.male = c(0)
fpr.test.male = c(0)

for(i in sort(test.result.p.male,decreasing = T)){
  tpr.test.male = c(tpr.test.male, sum(test.result.p.male >= i)/length(test.result.p.male))
  fpr.test.male = c(fpr.test.male, sum(test.result.h.male >= i)/length(test.result.h.male))
}

AUC.test.male = c()
for(i in 2:length(fpr.test.male)){
  S = (fpr.test.male[i]-fpr.test.male[i-1]) * (tpr.test.male[i-1] + (tpr.test.male[i]-tpr.test.male[i-1])/2 )
  AUC.test.male = sum(AUC.test.male, S)
}

################################ female
data.test.female <- data.test[data.test$sex=="女", c("age", "Alb", "CEA", "HCT", "LYMPH%", "class")]
pred.test.female <- predict(fit.CART.zl, data.test.female[, -6], type="prob")

table(data.test.female$class)
test.result.p.female = pred.test.female[1:table(data.test.female$class)[2],2]
test.result.h.female = pred.test.female[(table(data.test.female$class)[2]+1):dim(data.test.female)[1], 2]

tpr.test.female = c(0)
fpr.test.female = c(0)

for(i in sort(test.result.p.female,decreasing = T)){
  tpr.test.female = c(tpr.test.female, sum(test.result.p.female >= i)/length(test.result.p.female))
  fpr.test.female = c(fpr.test.female, sum(test.result.h.female >= i)/length(test.result.h.female))
}

AUC.test.female = c()
for(i in 2:length(fpr.test.female)){
  S = (fpr.test.female[i]-fpr.test.female[i-1]) * (tpr.test.female[i-1] + (tpr.test.female[i]-tpr.test.female[i-1])/2 )
  AUC.test.female = sum(AUC.test.female, S)
}


#####

#############################################################
##### (11) CART test (ratio)
##########################################################
############### rate01 842 : 10000
set.seed(4000)
data.label.rate.h = sample(seq(1,dim.test.h),10000)
data.test.h10000 = data.test.h[data.label.rate.h, ]
data.test.rate01 = rbind(data.test.p, data.test.h10000)
data.test.rate01 = data.test.rate01[, c(2,4,9,13,18,41)]

pred.CART.rate01 <- predict(fit.CART.zl, data.test.rate01, type="prob")

table(data.test.rate01$class)
test.result.p.rate01 = pred.CART.rate01[1:table(data.test.rate01$class)[2],2]
test.result.h.rate01 = pred.CART.rate01[(table(data.test.rate01$class)[2]+1):dim(data.test.rate01)[1], 2]

tpr.test.rate01 = c(0)
fpr.test.rate01 = c(0)

for(i in sort(test.result.p.rate01,decreasing = T)){
  tpr.test.rate01 = c(tpr.test.rate01, sum(test.result.p.rate01 >= i)/length(test.result.p.rate01))
  fpr.test.rate01 = c(fpr.test.rate01, sum(test.result.h.rate01 >= i)/length(test.result.h.rate01))
}

AUC.test.rate01 = c()
for(i in 2:length(fpr.test.rate01)){
  S = (fpr.test.rate01[i]-fpr.test.rate01[i-1]) * (tpr.test.rate01[i-1] + (tpr.test.rate01[i]-tpr.test.rate01[i-1])/2 )
  AUC.test.rate01 = sum(AUC.test.rate01, S)
}

sort(fpr.test.rate01, decreasing=FALSE)[471:480] ## 0.0081 0.0102
x1 = 0.0081; x2 = 0.0102
sort(tpr.test.rate01, decreasing=FALSE)[471:480] ## 0.5676960 0.6211401
y1 = 0.5676960; y2 = 0.6211401

tpr.rate01 = (0.01-x1)/(x2-x1)*(y2-y1)+y1




############### rate02 500 : 15419
set.seed(5000)
data.label.rate.p = sample(seq(1,dim.test.p),500)
data.test.p500 = data.test.p[data.label.rate.p, ]
data.test.rate02 = rbind(data.test.p500, data.test.h)
data.test.rate02 = data.test.rate02[, c(2,4,9,13,18,41)]

pred.CART.rate02 <- predict(fit.CART.zl, data.test.rate02, type="prob")

table(data.test.rate02$class)
test.result.p.rate02 = pred.CART.rate02[1:table(data.test.rate02$class)[2],2]
test.result.h.rate02 = pred.CART.rate02[(table(data.test.rate02$class)[2]+1):dim(data.test.rate02)[1], 2]

tpr.test.rate02 = c(0)
fpr.test.rate02 = c(0)

for(i in sort(test.result.p.rate02,decreasing = T)){
  tpr.test.rate02 = c(tpr.test.rate02, sum(test.result.p.rate02 >= i)/length(test.result.p.rate02))
  fpr.test.rate02 = c(fpr.test.rate02, sum(test.result.h.rate02 >= i)/length(test.result.h.rate02))
}

AUC.test.rate02 = c()
for(i in 2:length(fpr.test.rate02)){
  S = (fpr.test.rate02[i]-fpr.test.rate02[i-1]) * (tpr.test.rate02[i-1] + (tpr.test.rate02[i]-tpr.test.rate02[i-1])/2 )
  AUC.test.rate02 = sum(AUC.test.rate02, S)
}


sort(fpr.test.rate02, decreasing=FALSE)[326:330] ## 0.008949997 0.010441663
x1 = 0.008949997; x2 = 0.010441663
sort(tpr.test.rate02, decreasing=FALSE)[326:330] ## 0.656 0.710
y1 = 0.656; y2 = 0.710
tpr.rate02 = (0.01-x1)/(x2-x1)*(y2-y1)+y1

############### rate03 300 : 15419
set.seed(6000)
data.label.rate.p = sample(seq(1,dim.test.p),300)
data.test.p300 = data.test.p[data.label.rate.p, ]
data.test.rate03 = rbind(data.test.p300, data.test.h)
data.test.rate03 = data.test.rate03[, c(2,4,9,13,18,41)]

pred.CART.rate03 <- predict(fit.CART.zl, data.test.rate03, type="prob")

table(data.test.rate03$class)
test.result.p.rate03 = pred.CART.rate03[1:table(data.test.rate03$class)[2],2]
test.result.h.rate03 = pred.CART.rate03[(table(data.test.rate03$class)[2]+1):dim(data.test.rate03)[1], 2]

tpr.test.rate03 = c(0)
fpr.test.rate03 = c(0)

for(i in sort(test.result.p.rate03,decreasing = T)){
  tpr.test.rate03 = c(tpr.test.rate03, sum(test.result.p.rate03 >= i)/length(test.result.p.rate03))
  fpr.test.rate03 = c(fpr.test.rate03, sum(test.result.h.rate03 >= i)/length(test.result.h.rate03))
}

AUC.test.rate03 = c()
for(i in 2:length(fpr.test.rate03)){
  S = (fpr.test.rate03[i]-fpr.test.rate03[i-1]) * (tpr.test.rate03[i-1] + (tpr.test.rate03[i]-tpr.test.rate03[i-1])/2 )
  AUC.test.rate03 = sum(AUC.test.rate03, S)
}


sort(fpr.test.rate03, decreasing=FALSE)[186:195] ## 0.008949997 0.010441663
x1 = 0.008949997; x2 = 0.010441663
sort(tpr.test.rate03, decreasing=FALSE)[186:195] ## 0.63 0.69
y1 = 0.63; y2 = 0.69
tpr.rate03 = (0.01-x1)/(x2-x1)*(y2-y1)+y1
tpr.rate03

############### rate04 100 : 15419
set.seed(7000)
data.label.rate.p = sample(seq(1,dim.test.p),100)
data.test.p100 = data.test.p[data.label.rate.p, ]
data.test.rate04 = rbind(data.test.p100, data.test.h)
data.test.rate04 = data.test.rate04[, c(2,4,9,13,18,41)]

pred.CART.rate04 <- predict(fit.CART.zl, data.test.rate04, type="prob")

table(data.test.rate04$class)
test.result.p.rate04 = pred.CART.rate04[1:table(data.test.rate04$class)[2],2]
test.result.h.rate04 = pred.CART.rate04[(table(data.test.rate04$class)[2]+1):dim(data.test.rate04)[1], 2]

tpr.test.rate04 = c(0)
fpr.test.rate04 = c(0)

for(i in sort(test.result.p.rate04,decreasing = T)){
  tpr.test.rate04 = c(tpr.test.rate04, sum(test.result.p.rate04 >= i)/length(test.result.p.rate04))
  fpr.test.rate04 = c(fpr.test.rate04, sum(test.result.h.rate04 >= i)/length(test.result.h.rate04))
}

AUC.test.rate04 = c()
for(i in 2:length(fpr.test.rate04)){
  S = (fpr.test.rate04[i]-fpr.test.rate04[i-1]) * (tpr.test.rate04[i-1] + (tpr.test.rate04[i]-tpr.test.rate04[i-1])/2 )
  AUC.test.rate04 = sum(AUC.test.rate04, S)
}

sort(fpr.test.rate04, decreasing=FALSE)[61:65] ## 0.008949997 0.010441663
x1 = 0.008949997; x2 =0.010441663
sort(tpr.test.rate04, decreasing=FALSE)[61:65] ## 0.63 0.69
y1 = 0.63; y2 = 0.69
tpr.rate04 = (0.01-x1)/(x2-x1)*(y2-y1)+y1
tpr.rate04


############### rate05 50 : 15419
set.seed(8000)
data.label.rate.p = sample(seq(1,dim.test.p),50)
data.test.p50 = data.test.p[data.label.rate.p, ]
data.test.rate05 = rbind(data.test.p50, data.test.h)
data.test.rate05 = data.test.rate05[, c(2,4,9,13,18,41)]

pred.CART.rate05 <- predict(fit.CART.zl, data.test.rate05, type="prob")

table(data.test.rate05$class)
test.result.p.rate05 = pred.CART.rate05[1:table(data.test.rate05$class)[2],2]
test.result.h.rate05 = pred.CART.rate05[(table(data.test.rate05$class)[2]+1):dim(data.test.rate05)[1], 2]

tpr.test.rate05 = c(0)
fpr.test.rate05 = c(0)

for(i in sort(test.result.p.rate05,decreasing = T)){
  tpr.test.rate05 = c(tpr.test.rate05, sum(test.result.p.rate05 >= i)/length(test.result.p.rate05))
  fpr.test.rate05 = c(fpr.test.rate05, sum(test.result.h.rate05 >= i)/length(test.result.h.rate05))
}

AUC.test.rate05 = c()
for(i in 2:length(fpr.test.rate05)){
  S = (fpr.test.rate05[i]-fpr.test.rate05[i-1]) * (tpr.test.rate05[i-1] + (tpr.test.rate05[i]-tpr.test.rate05[i-1])/2 )
  AUC.test.rate05 = sum(AUC.test.rate05, S)
}


sort(fpr.test.rate05, decreasing=FALSE)[31:35] ## 0.008949997 0.010441663
x1 = 0.008949997; x2 = 0.010441663
sort(tpr.test.rate05, decreasing=FALSE)[31:35] ## 0.63 0.69
y1 = 0.62; y2 = 0.70
tpr.rate05 = (0.01-x1)/(x2-x1)*(y2-y1)+y1
tpr.rate05


############### rate06 30 : 15419
set.seed(9000)
data.label.rate.p = sample(seq(1,dim.test.p),30)
data.test.p30 = data.test.p[data.label.rate.p, ]
data.test.rate06 = rbind(data.test.p30, data.test.h)
data.test.rate06 = data.test.rate06[, c(2,4,9,13,18,41)]

pred.CART.rate06 <- predict(fit.CART.zl, data.test.rate06, type="prob")

table(data.test.rate06$class)
test.result.p.rate06 = pred.CART.rate06[1:table(data.test.rate06$class)[2],2]
test.result.h.rate06 = pred.CART.rate06[(table(data.test.rate06$class)[2]+1):dim(data.test.rate06)[1], 2]

tpr.test.rate06 = c(0)
fpr.test.rate06 = c(0)

for(i in sort(test.result.p.rate06,decreasing = T)){
  tpr.test.rate06 = c(tpr.test.rate06, sum(test.result.p.rate06 >= i)/length(test.result.p.rate06))
  fpr.test.rate06 = c(fpr.test.rate06, sum(test.result.h.rate06 >= i)/length(test.result.h.rate06))
}

AUC.test.rate06 = c()
for(i in 2:length(fpr.test.rate06)){
  S = (fpr.test.rate06[i]-fpr.test.rate06[i-1]) * (tpr.test.rate06[i-1] + (tpr.test.rate06[i]-tpr.test.rate06[i-1])/2 )
  AUC.test.rate06 = sum(AUC.test.rate06, S)
}


sort(fpr.test.rate06, decreasing=FALSE) ## 0.008949997 0.010441663
x1 = 0.008949997; x2 = 0.010441663
sort(tpr.test.rate06, decreasing=FALSE) ## 0.63 0.69
y1 = 0.6666667; y2 = 0.70
tpr.rate06 = (0.01-x1)/(x2-x1)*(y2-y1)+y1
tpr.rate06


############### rate07 10 : 15419
set.seed(10000)
data.label.rate.p = sample(seq(1,dim.test.p),10)
data.test.p10 = data.test.p[data.label.rate.p, ]
data.test.rate07 = rbind(data.test.p10, data.test.h)
data.test.rate07 = data.test.rate07[, c(2,4,9,13,18,41)]

pred.CART.rate07 <- predict(fit.CART.zl, data.test.rate07, type="prob")

table(data.test.rate07$class)
test.result.p.rate07 = pred.CART.rate07[1:table(data.test.rate07$class)[2],2]
test.result.h.rate07 = pred.CART.rate07[(table(data.test.rate07$class)[2]+1):dim(data.test.rate07)[1], 2]

tpr.test.rate07 = c(0)
fpr.test.rate07 = c(0)

for(i in sort(test.result.p.rate07,decreasing = T)){
  tpr.test.rate07 = c(tpr.test.rate07, sum(test.result.p.rate07 >= i)/length(test.result.p.rate07))
  fpr.test.rate07 = c(fpr.test.rate07, sum(test.result.h.rate07 >= i)/length(test.result.h.rate07))
}

AUC.test.rate07 = c()
for(i in 2:length(fpr.test.rate07)){
  S = (fpr.test.rate07[i]-fpr.test.rate07[i-1]) * (tpr.test.rate07[i-1] + (tpr.test.rate07[i]-tpr.test.rate07[i-1])/2 )
  AUC.test.rate07 = sum(AUC.test.rate07, S)
}


sort(fpr.test.rate07, decreasing=FALSE) ## 0.008949997 0.010441663
x1 = 0.008949997; x2 = 0.010441663
sort(tpr.test.rate07, decreasing=FALSE) ## 0.63 0.69
y1 = 0.6666667; y2 = 0.70
tpr.rate07 = (0.01-x1)/(x2-x1)*(y2-y1)+y1
tpr.rate07 = 0.85

############### rate08 4 : 15419
set.seed(21000)
data.label.rate.p = sample(seq(1,dim.test.p), 4)
data.test.p4 = data.test.p[data.label.rate.p, ]
data.test.rate08 = rbind(data.test.p4, data.test.h)
data.test.rate08 = data.test.rate08[, c(2,4,9,13,18,41)]

pred.CART.rate08 <- predict(fit.CART.zl, data.test.rate08, type="prob")

table(data.test.rate08$class)
test.result.p.rate08 = pred.CART.rate08[1:table(data.test.rate08$class)[2],2]
test.result.h.rate08 = pred.CART.rate08[(table(data.test.rate08$class)[2]+1):dim(data.test.rate08)[1], 2]

tpr.test.rate08 = c(0)
fpr.test.rate08 = c(0)

for(i in sort(test.result.p.rate08,decreasing = T)){
  tpr.test.rate08 = c(tpr.test.rate08, sum(test.result.p.rate08 >= i)/length(test.result.p.rate08))
  fpr.test.rate08 = c(fpr.test.rate08, sum(test.result.h.rate08 >= i)/length(test.result.h.rate08))
}

AUC.test.rate08 = c()
for(i in 2:length(fpr.test.rate08)){
  S = (fpr.test.rate08[i]-fpr.test.rate08[i-1]) * (tpr.test.rate08[i-1] + (tpr.test.rate08[i]-tpr.test.rate08[i-1])/2 )
  AUC.test.rate08 = sum(AUC.test.rate08, S)
}


sort(fpr.test.rate08, decreasing=FALSE) ## 0.008949997 0.010441663
x1 = 0.007263765; x2 = 1
sort(tpr.test.rate08, decreasing=FALSE) ## 0.63 0.69
y1 = 0.75; y2 = 1
tpr.rate08 = (0.01-x1)/(x2-x1)*(y2-y1)+y1
tpr.rate08




#### 处理

AUC.rate = c(AUC.test.rate01, AUC.CART.zl.test, AUC.test.rate02, AUC.test.rate03, AUC.test.rate04,
             AUC.test.rate05, AUC.test.rate06, AUC.test.rate07, AUC.test.rate08
             )


tpr.rate = c(tpr.rate01, 0.73, tpr.rate02, tpr.rate03, tpr.rate04, tpr.rate05, tpr.rate06,
             tpr.rate07, tpr.rate08)



######

#############################################################
##### (12) sg
##########################################################
###  预处理
head(patient.info.sg)
colnames(patient.info.sg)
patient.matrix.sg = cbind( patient.info.sg[, c(5,24,22,27,20)], 
                           class = rep( "patient", dim(patient.info.sg)[1] ),
                           stringsAsFactors = FALSE
)
colnames(patient.matrix.sg) = c("age", "Alb", "CEA", "HCT", "LYMPH%", "class")
head(patient.matrix.sg)

head(health.info.sg)
colnames(health.info.sg)
health.matrix.sg  = cbind( health.info.sg[, c(3,22,20,25,18)], 
                           class = rep( "health", dim(health.info.sg)[1] ),
                           stringsAsFactors = FALSE
                           )
colnames(health.matrix.sg) = c("age", "Alb", "CEA", "HCT", "LYMPH%", "class")
head(health.matrix.sg)

colnames(patient.matrix.sg)
colnames(health.matrix.sg)

dim.sg.p = dim(patient.matrix.sg)[1]
dim.sg.h = dim(health.matrix.sg)[1]

dim.sgtest.p = floor(dim.sg.p*0.3)
dim.sgtest.h = floor(dim.sg.h*0.3)

set.seed(30000)
data.label.sgtest.p  = sample( seq(1,dim.sg.p), dim.sgtest.p, replace=FALSE )
data.label.sgtrain.p = setdiff( seq(1,dim.sg.p), data.label.sgtest.p )

set.seed(40000)
data.label.sgtest.h  = sample( seq(1,dim.sg.h), dim.sgtest.h, replace=FALSE )
data.label.sgtrain.h = setdiff( seq(1,dim.sg.h), data.label.sgtest.h )

data.sgtest.p = patient.matrix.sg[data.label.sgtest.p, ]
data.sgtest.h = health.matrix.sg[data.label.sgtest.h, ]
data.sgtest   = rbind( data.sgtest.p, data.sgtest.h )

data.sgtrain.p = patient.matrix.sg[data.label.sgtrain.p, ]
data.sgtrain.h = health.matrix.sg[data.label.sgtrain.h, ]
data.sgtrain   = rbind( data.sgtrain.p, data.sgtrain.h )

### fit.CART.zl, sgtest
pred.CART.zl.sg <- predict(fit.CART.zl, data.sgtest[,-6], type="prob")

table(data.sgtest$class)
test.result.p.zl.sg = pred.CART.zl.sg[1:table(data.sgtest$class)[2],2]
test.result.h.zl.sg = pred.CART.zl.sg[(table(data.sgtest$class)[2]+1):dim(data.sgtest)[1], 2]

tpr.test.zl.sg = c(0)
fpr.test.zl.sg = c(0)
test.result.zl.sg = unique(c(test.result.p.zl.sg,test.result.h.zl.sg))
for(i in sort(test.result.zl.sg,decreasing=TRUE)){
  tpr.test.zl.sg = c(tpr.test.zl.sg, sum(test.result.p.zl.sg >= i)/length(test.result.p.zl.sg))
  fpr.test.zl.sg = c(fpr.test.zl.sg, sum(test.result.h.zl.sg >= i)/length(test.result.h.zl.sg))
}

AUC.test.zl.sg = c()
for(i in 2:length(fpr.test.zl.sg)){
  S = (fpr.test.zl.sg[i]-fpr.test.zl.sg[i-1]) * (tpr.test.zl.sg[i-1] + (tpr.test.zl.sg[i]-tpr.test.zl.sg[i-1])/2 )
  AUC.test.zl.sg = sum(AUC.test.zl.sg, S)
}


### fit.CART.sg, sgtest
fit.CART.sg = rpart(class~., data=data.sgtrain)

pred.CART.sg <- predict(fit.CART.sg, data.sgtest[,-6], type="prob")

table(data.sgtest$class)
test.result.p.sg = pred.CART.sg[1:table(data.sgtest$class)[2],2]
test.result.h.sg = pred.CART.sg[(table(data.sgtest$class)[2]+1):dim(data.sgtest)[1], 2]

tpr.test.sg = c(0)
fpr.test.sg = c(0)
test.result.sg = unique(c(test.result.p.sg,test.result.h.sg))
for(i in sort(test.result.sg,decreasing=TRUE)){
  tpr.test.sg = c(tpr.test.sg, sum(test.result.p.sg >= i)/length(test.result.p.sg))
  fpr.test.sg = c(fpr.test.sg, sum(test.result.h.sg >= i)/length(test.result.h.sg))
}

AUC.test.sg = c()
for(i in 2:length(fpr.test.sg)){
  S = (fpr.test.sg[i]-fpr.test.sg[i-1]) * (tpr.test.sg[i-1] + (tpr.test.sg[i]-tpr.test.sg[i-1])/2 )
  AUC.test.sg = sum(AUC.test.sg, S)
}


#####



## ##########################################################
save.image(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/analysis/zl_hospital_analysis.RData"
  )

#####################################################################
##### (7)随机打乱标签
##########################################################

###  保存好原数据
tmp.p = patient.matrix.zl
tmp.h = health.matrix.zl

## 预处理
all.matrix = rbind(patient.matrix.zl, health.matrix.zl)
set.seed(1)
patient.label = sample( 1:dim(all.matrix)[1], 4211 )
health.label  = setdiff( seq(1,dim(all.matrix)[1]), patient.label )

patient.matrix.zl = all.matrix[patient.label, ]
health.matrix.zl  = all.matrix[health.label, ]

patient.matrix.zl$class = "patient"
health.matrix.zl$class  = "health"

##  重复上面的过程(4)

fit.permutate = rpart(class~., data=data.train)

tpr.CART.permute = c(0)
fpr.CART.permute = c(0)

for(i in sort(valid.result.p,decreasing = T)){
  tpr.CART.permute = c(tpr.CART.permute, sum(valid.result.p >= i)/length(valid.result.p))
  fpr.CART.permute = c(fpr.CART.permute, sum(valid.result.h >= i)/length(valid.result.h))
}


##  复原数据
patient.matrix.zl = tmp.p
health.matrix.zl  = tmp.h

rm(tmp.p,tmp.h,all.matrix)

#####

