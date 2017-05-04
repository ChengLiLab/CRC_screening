# data preprocess of the health data

setwd("/lustre/user/liclab/shiq/CRC/sg_hospital/health")

###################################################################################
### (1)import datasets
# library( data.table )
health.rawdata.sg <- read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/health/20160729前体检数据.csv",
  header = T, stringsAsFactors = F )

dim( health.rawdata.sg ) ## 5709152      15
length( unique( health.rawdata.sg$MEDICAL_NO ) ) ## 66570 number of samples
length( unique( health.rawdata.sg$ITEM_NAME ) )  ## 434 items
range( table( health.rawdata.sg$MEDICAL_NO ) )  ## 1 688 range of items


##################################################################################
################################ (2) extract the useful columns
health.data.sg = health.rawdata.sg[, c("MEDICAL_NO", "SEX", "AGE", "ITEM_NAME", "MEDICAL_VALUE", "REFERENCE_RANGE") ]
colnames( health.data.sg ) = c( "no", "sex", "age", "item.c", "value", "ref" )
health.data.sg = unique( health.data.sg )
dim( health.data.sg ) ## 5647709       6

length( unique( health.data.sg$no ) ) ## 66570 number of samples
length( unique( health.data.sg$item.c ) ) ## 434 items
range( table( health.data.sg$no ) )  ## 1 207 range of items

##############################################################################
################## (3)examine each columns
head(table( health.data.sg$no )) ## 无no为空
tail(table( health.data.sg$no ))
sum( table(health.data.sg$no) )
sum( health.data.sg$no == "" )  ## 0
length( unique( health.data.sg$no) ) ## 66570 number of samples

table( health.data.sg$sex )
sum( table(health.data.sg$sex) )

table( health.data.sg$age )
## 查看age = 2015的体检日期,2015/7/17"
health.rawdata.sg[health.rawdata.sg$AGE == "2015", "MEDICAL_DATE"] 
## 查看age = 2015的生日,0076/3/3
health.rawdata.sg[health.rawdata.sg$AGE == "2015", "BIRTHDAY"]    
health.data.sg$age[health.data.sg$age == 2015 ] = 39

table( health.data.sg$item.c )
unique( health.data.sg$item.c ) ## 434 items
table(health.data.sg$item.c[grep("百分比", health.data.sg$item.c)])
health.data.sg$item.c[health.data.sg$item.c == "单核细胞百分比"] = "单核细胞比率"

health.data.sg[(health.data.sg$item.c == "白细胞") & (health.data.sg$ref == "阴性"),
               "item.c"] = "尿白细胞"
health.data.sg$item.c[health.data.sg$item.c == "白细胞数"] = "白细胞"

health.data.sg$item.c[ grep( "^[ ]", health.data.sg$item.c ) ] =
  gsub( "^[ ]", "", health.data.sg$item.c[grep("^[ ]", health.data.sg$item.c)] )

head(table( health.data.sg$value ))
tail(table( health.data.sg$value ))
sum( table( health.data.sg$value ) ) ## 5647709
sum( is.na( health.data.sg$value ) ) ## 0

health.data.sg$value[health.data.sg$value == 1000]
health.data.sg$value[health.data.sg$value == 2000]

##  先把阴阳保护起来，可以继续推敲
health.data.sg$value[ health.data.sg$value == "阴性" ]  = 1000
health.data.sg$value[ health.data.sg$value == "-----" ] = 1000
health.data.sg$value[ health.data.sg$value == "----" ]  = 1000
health.data.sg$value[ health.data.sg$value == "--" ]    = 1000
health.data.sg$value[ health.data.sg$value == "-" ]     = 1000

health.data.sg$value[ health.data.sg$value == "阳性" ]  = 2000
health.data.sg$value[ health.data.sg$value == "+" ]     = 2000
health.data.sg$value[ health.data.sg$value == "++" ]    = 2000
health.data.sg$value[ health.data.sg$value == "+++" ]   = 2000

health.data.sg$value = as.numeric( health.data.sg$value )

health.data.sg$value[ health.data.sg$value == 1000 ] = "阴"
health.data.sg$value[ health.data.sg$value == 2000 ] = "阳"
## 问题：数值型指标的value有的是用阴阳标记的，这部分也被保护了
## 解决：最后将将数值型指标再as.numeric(),阴阳变成了NA 

table( health.data.sg$ref )

### delete the data,whose "no","sex","age",or "value" are empty
health.data.sg = na.omit( health.data.sg, cols = c("no", "sex", "age", "item", "value") )

head(health.data.sg)
dim( health.data.sg ) ## 3247030       6
length( unique( health.data.sg$no) ) ## 63520 number of samples
length( unique( health.data.sg$item.c) )  ## 290 items

### delete repeating data,just reserve the newest one
head(health.data.sg)
health.data.sg.paste = paste0( health.data.sg$no, "_", 
                            health.data.sg$item )
id          = match( unique(health.data.sg.paste), health.data.sg.paste )
health.data.sg = health.data.sg[ id, ]

rownames(health.data.sg) = NULL

length( unique( health.data.sg$no ) ) ## 63520 number of samples
length( unique( health.data.sg$item.c ) ) ## 290
range( table( health.data.sg$no ) )  ## 1 158 range of items
unique( health.data.sg$item ) ## 290 items

write.table( health.data.sg,
             "/lustre/user/liclab/shiq/CRC/sg_hospital/health/health.data.sg.txt",
             quote = T )

## ############################################################################ 
####### (4)结合patient与health数据确定研究指标
length( table(health.data.sg$item) ) ## 290 items

sort( table( health.data.sg$item.c ),decreasing=TRUE)
# sum( table(health.data.sg$item) >= 15000 ) ## 80 items whose samples are more than 15000

health.item.40.sg = rownames( 
  sort(
    table(health.data.sg$item),decreasing=TRUE
  )
)[1: floor((length(unique(health.data.sg$item)))*0.4)] ## 59


write.table( health.item.40.sg,
             "/lustre/user/liclab/shiq/CRC/sg_hospital/health/health.item.40.sg.txt",
             quote = T )

patient.item.35.sg = read.table(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/patient.item.35.sg.txt",
  stringsAsFactors = F
  ) ## 71 items

patient.item.35.sg = c(patient.item.35.sg$x)

items.sg = intersect( patient.item.35.sg, health.item.40.sg ) ## 28 items

write.table( items.sg,
             "/lustre/user/liclab/shiq/CRC/sg_hospital/health/items.sg.txt",
             quote = T )

unique(sort(table(health.data.sg$item)))

############################################################################
## (5) 筛选 28 items 的数据 
health.data.sg.new = as.data.frame(matrix(numeric(0),ncol=4)) ## data about 46 items
health.ref.sg = as.data.frame(matrix(numeric(0),ncol=2))
for(i in items.sg){
  health.data.sg.new.single = health.data.sg[ health.data.sg$item == i, ]
  health.ref.sg             = rbind(health.ref.sg, health.data.sg.new.single[1,c("item.c","ref")])
  health.data.sg.new        = rbind( health.data.sg.new, health.data.sg.new.single )
}

health.data.sg.new = health.data.sg.new[order(health.data.sg.new$no), ]
rownames(health.data.sg.new) = NULL
rownames(health.ref.sg)      = NULL
rm( health.data.sg.new.single )

length( unique( health.data.sg.new$no ) ) ## 56799 number of samples
unique( health.data.sg.new$item ) ## 28 items

write.table(health.ref.sg, 
            "/lustre/user/liclab/shiq/CRC/sg_hospital/health/health.ref.sg.txt",
            quote = T, sep = "\t"
)

write.table(health.data.sg.new, 
            "/lustre/user/liclab/shiq/CRC/sg_hospital/health/health.data.sg.new.txt",
            quote = F
)


############################################################################
####### transform raw data "health.data.sg_new" into "health.matrix"

library(reshape2)
health.info.sg = dcast(health.data.sg.new, no+sex+age  ~ item.c, value = "value")

###  查看指标是定性还是定量

#### 将变量值由character变成numeric
items.num = setdiff( colnames(health.info.sg), c("no", "sex") )
for (i in items.num) {
  health.info.sg[, i] = as.numeric( health.info.sg[, i] )
}

write.table(health.info.sg, 
            "/lustre/user/liclab/shiq/CRC/sg_hospital/health/health.info.sg.txt", 
            quote = T
)

############################################################################

save.image(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/health/sg_hospital_health.RData"
)

