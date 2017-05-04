## preprocess of health data of zl hospital
## 2016.04, shiq

setwd("/lustre/user/liclab/shiq/CRC/zl_hospital/health")

## #################################################################################
####### (1) import data 2007-2014
health.2007 = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/07-14年含年龄性别体检/2007体检.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
health.2008 = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/07-14年含年龄性别体检/2008体检.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
health.2009 = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/07-14年含年龄性别体检/2009体检.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
health.2010 = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/07-14年含年龄性别体检/2010体检.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
health.2011 = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/07-14年含年龄性别体检/2011体检.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
health.2012 = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/07-14年含年龄性别体检/2012体检.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
health.2013 = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/07-14年含年龄性别体检/2013体检.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
health.2014 = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/07-14年含年龄性别体检/2014体检.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )

health.rawdata.zl = rbind( health.2007, 
                        health.2008,
                        health.2009,
                        health.2010,
                        health.2011,
                        health.2012, 
                        health.2013,
                        health.2014
                        )
dim( health.rawdata.zl ) ## 2442954       8
length( unique( health.rawdata.zl$p_no ) ) ## 80194 samples
length( unique( health.rawdata.zl$item ) ) ## 220 items

rm( health.2007,
    health.2008, 
    health.2009, 
    health.2010,
    health.2011, 
    health.2012, 
    health.2013, 
    health.2014
    )

## #################################################################################
##### (2) extract the useful columns

head(health.rawdata.zl)
colnames(health.rawdata.zl) = c( "no", "sex", "age", "item", "item.c", "value", "unit", "ref" )

health.data.zl = health.rawdata.zl[, c( "no", "sex", "age", "item", "value", "ref" ) ]
head( health.data.zl )

## ##############################################################################
### (3) examine each columns

head(table( health.data.zl$no )) ## 有no为空
tail(table( health.data.zl$no ))
sum( table(health.data.zl$no ))
length( unique( health.data.zl$no) ) ## 80194 number of samples
health.data.zl$no[ health.data.zl$no == "" ]     = NA # 119380
health.data.zl$no[ health.data.zl$no == "无" ]   = NA # 10
health.data.zl$no[ health.data.zl$no == "无号" ] = NA # 8
## delete "." of tail of some no.
health.data.zl$no[ grep( "[.]$", health.data.zl$no ) ] = 
  gsub( "[.]$", "", health.data.zl$no[ grep("[.]$", health.data.zl$no) ] )
length( unique( health.data.zl$no) ) ## 78744 number of samples


table( health.data.zl$sex )
sum( table( health.data.zl$sex ) )
health.data.zl$sex[ health.data.zl$sex == "南" ] = "男"
health.data.zl$sex[ health.data.zl$sex == "" ]   =  NA


table( health.data.zl$age )
sum( table( health.data.zl$age ) )
sum( is.na(health.data.zl$age)) ## 160619
## 将age > 110的视为录入错误，删除 
health.data.zl$age[ health.data.zl$age > 110 ] = NA ## 90
sum( is.na(health.data.zl$age) ) ## 160709


sum( table( health.data.zl$item ) )
unique( health.data.zl$item ) ## 220 items
sum( health.data.zl$item == "" ) ## 2
health.data.zl$item[ health.data.zl$item == "" ] = NA

## delete "." of tail of some items
# method 1
health.data.zl$item[ grep( "[.]$", health.data.zl$item ) ] = 
  gsub( "[.]$", "", health.data.zl$item[grep("[.]$", health.data.zl$item)] )
unique( health.data.zl$item ) ## 188 items
# method 2
# health.data.zl$item[ grep("[.]$", health.data.zl$item) ] = 
# substr(health.data.zl$item[grep( "[.]$", health.data.zl$item)], 
# start = c(1,1), 
# stop = nchar(health.data.zl$item[grep(pattern = "[.]$", x = health.data.zl$item)]) - 1)

## divide "RDW" into "RDW-SD" and "RDW-CV" by "ref"
health.data.zl$item[ (health.data.zl$item == "RDW") & (health.data.zl$ref == "37.0-50.0") ] =
  "RDW-SD"
health.data.zl$item[ (health.data.zl$item == "RDW") & 
                    ( health.data.zl$ref == "11.6-14.8" |  health.data.zl$ref == "11.60-14.80" )
                  ] =
  "RDW-CV"
health.data.zl$item[ health.data.zl$item == "T-Bil" ] = "TBil"

unique( health.data.zl$item ) ## 186 items


head(table( health.data.zl$value ))
tail(table( health.data.zl$value ))
sum(table( health.data.zl$value ))
sum( health.data.zl$value == "" )
##  先把阴阳保护起来, 可以继续推敲
health.data.zl$value[ grep( "[－]", health.data.zl$value ) ] = 1000
health.data.zl$value[ grep( "[—]", health.data.zl$value ) ]  = 1000
health.data.zl$value[ grep( "[_]", health.data.zl$value ) ]  = 1000

health.data.zl$value[ grep( "[阴]", health.data.zl$value ) ] = 1000
health.data.zl$value[ health.data.zl$value == "-性" ]        = 1000
health.data.zl$value[ health.data.zl$value == "NEGATIVE" ]   = 1000

health.data.zl$value[ grep( "[阳]", health.data.zl$value ) ] = 2000
health.data.zl$value[ grep( "[+]", health.data.zl$value ) ]  = 2000
health.data.zl$value[ grep( "[＋]", health.data.zl$value ) ] = 2000
health.data.zl$value[ health.data.zl$value == "POSITIVE" ]   = 2000

health.data.zl$value = as.numeric( health.data.zl$value )

health.data.zl$value[ health.data.zl$value == 1000 ] = "阴"
health.data.zl$value[ health.data.zl$value == 2000 ] = "阳"


table( health.data.zl$ref )
health.data.zl$ref[ grep( "[阴]", health.data.zl$ref ) ] = "阴性"

### delete the data,whose "no", "item", or "value" are empty
health.data.zl = na.omit( health.data.zl, cols = c("no", "sex", "age", "item", "value") )

dim( health.data.zl ) ## 2047406       6
length( unique( health.data.zl$no ) ) ## 77169 samples
length( unique( health.data.zl$item ) ) ## 169 items


### delete repeating data,just reserve the newest one
head(health.data.zl)
health.data.zl.paste = paste0( health.data.zl$no, "_", 
                            health.data.zl$item 
                            )
id          = match( unique(health.data.zl.paste), health.data.zl.paste )
health.data.zl = health.data.zl[ id, ]

rownames(health.data.zl) = NULL

dim( health.data.zl ) ## 1325486       6
length( unique( health.data.zl$no ) ) ## 77169 number of samples
range( table( health.data.zl$no ) )  ## 1 97 range of items
unique( health.data.zl$item ) ## 169 items

write.table( health.data.zl,
             "/lustre/user/liclab/shiq/CRC/zl_hospital/health/health.data.zl.txt",
             quote = T )

## #########################################################################
####### (4)结合patient与health数据确定研究指标
length( unique( health.data.zl$no ) ) ## 77169 number of samples
range( table( health.data.zl$no ) )  ## 1 97 range of items
length( table(health.data.zl$item) ) ## 169 items

sort( table( health.data.zl$item ), decreasing=TRUE )
# sum( table(health.data.zl$item) >= 1400 ) ## 67 items whose samples are more than 1400

health.item.40.zl = rownames( 
  sort(
    table(health.data.zl$item),decreasing=TRUE
  )
)[1: floor((length(unique(health.data.zl$item)))*0.4)] ## 59

write.table( health.item.40.zl,
             "/lustre/user/liclab/shiq/CRC/zl_hospital/health/health.item.40.zl.txt",
             quote = T )

patient.item.35.zl = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/patient.item.35.zl.txt",
  stringsAsFactors = F
  ) ## 87 items

patient.item.35.zl = c(patient.item.35.zl$x)

items.zl = intersect( patient.item.35.zl, health.item.40.zl ) ## 38 items

write.table( items.zl,
             "/lustre/user/liclab/shiq/CRC/zl_hospital/health/items.zl.txt",
             quote = T )


############################################################################
### (5) 筛选 38 items 的数据
health.data.zl.new = as.data.frame(matrix(numeric(0),ncol=4)) ## data about 38 items
health.ref.zl = as.data.frame(matrix(numeric(0),ncol=2))
for(i in 1:length(items.zl)){
  health.data.zl.new.single = health.data.zl[ health.data.zl$item == items.zl[i], ]
  health.ref.zl             = rbind(health.ref.zl, health.data.zl.new.single[1,c("item","ref")])
  health.data.zl.new        = rbind( health.data.zl.new, health.data.zl.new.single )
  }

health.data.zl.new = health.data.zl.new[order(health.data.zl.new$no), ]
rownames(health.data.zl.new) = NULL
rownames(health.ref.zl)      = NULL
rm( health.data.zl.new.single )

dim(health.data.zl.new) ## 1013651       6
length( unique( health.data.zl.new$no ) ) ## 77099 number of samples
range( table( health.data.zl.new$no ) )  ## 1 38 range of items
unique( health.data.zl.new$item ) ## 40 items

write.table(health.ref.zl, 
            "/lustre/user/liclab/shiq/CRC/zl_hospital/health/health.ref.zl.txt",
            quote = T, sep = "\t"
            )

write.table(health.data.zl.new, 
            "/lustre/user/liclab/shiq/CRC/zl_hospital/health/health.data.zl.new.txt",
            quote = F
            )

############################################################################
####### transform raw data "health.data.zl_new" into "health.matrix"
library(reshape2)
tmp = dcast( health.data.zl.new[, -c(2,3)], no ~ item, value = "value" )

####### correct the wrong age and sex
id  = match( tmp$no, health.data.zl.new$no )
sex = health.data.zl.new$sex[id]
age = health.data.zl.new$age[id]

head(tmp)
health.info.zl = cbind( no = tmp[, "no"], 
                        sex,
                        age,
                        tmp[, -c(1)], 
                        stringsAsFactors = FALSE
                        )
head(health.info.zl)

###  查看指标是定性还是定量
###  HBeAg HBsAg在health中基本全是定性，而在patient中是定量，且ref不一致，故删除
health.info.zl = health.info.zl[, setdiff( colnames(health.info.zl), c("HBeAg", "HBsAg")) ]

### as.numeric
items.num = setdiff( colnames(health.info.zl), c("no", "sex") )
for (i in items.num) {
  health.info.zl[, i] = as.numeric( health.info.zl[, i] )
}

write.table(health.info.zl, 
            "/lustre/user/liclab/shiq/CRC/zl_hospital/health/health.info.zl.txt", 
            quote = T
            )
rm(tmp)

############################################################################
save.image(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/zl_hospital_health.RData"
  )

