# preprocess of patient data of sg hospital 2011-20160826
## 2016.04, shiq

setwd("/lustre/user/liclab/shiq/CRC/sg_hospital/patient")

## #################################################################################
####### (1) import data 011-20160826
patient.colon.2011_2013 = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/20160826结直肠/结肠癌18/结肠11-13.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
patient.colon.2014 = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/20160826结直肠/结肠癌18/14.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.colon.2015 = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/20160826结直肠/结肠癌18/15.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.colon.2016 = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/20160826结直肠/结肠癌18/16.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)

patient.rectum.2011_2013 = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/20160826结直肠/直肠癌20/直肠11-13.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.rectum.2014 = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/20160826结直肠/直肠癌20/14.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.rectum.2015 = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/20160826结直肠/直肠癌20/15.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.rectum.2016 = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/20160826结直肠/直肠癌20/16.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)


head(patient.colon.2011_2013)
head(patient.colon.2014)
head(patient.colon.2015)

### 删除最后一列
patient.colon.2011_2013  = patient.colon.2011_2013[, -18]
patient.colon.2014       = patient.colon.2014[, -18]
patient.rectum.2011_2013 = patient.rectum.2011_2013[, -18]

### 重新统一命名
col.names = c( "no", "name","age", "sex", 
               "start_date_time","admission_date_time", "discharge_date_time",
               "stage_T","stage_N","stage_M","stages",
               "specimen",
               "result_date_time", "item.c", "value", "unit","ref"
               )
colnames(patient.colon.2011_2013)  = col.names
colnames(patient.colon.2014)       = col.names
colnames(patient.colon.2015)       = col.names
colnames(patient.colon.2016)       = col.names

colnames(patient.rectum.2011_2013) = col.names
colnames(patient.rectum.2014)      = col.names
colnames(patient.rectum.2015)      = col.names
colnames(patient.rectum.2016)      = col.names

### 添加位置信息
patient.colon  = rbind( patient.colon.2011_2013,
                        patient.colon.2014,
                        patient.colon.2015,
                        patient.colon.2016
                        )
patient.rectum = rbind( patient.rectum.2011_2013,
                        patient.rectum.2014,
                        patient.rectum.2015,
                        patient.rectum.2016
                        )

patient.colon  = cbind( patient.colon[, c("no","name")],
                        location = rep( "colon", dim(patient.colon)[1] ),
                        patient.colon[, -c(1,2)]
                        )
patient.rectum = cbind( patient.rectum[, c("no","name")],
                        location = rep( "rectum", dim(patient.rectum)[1] ),
                        patient.rectum[, -c(1,2)]
                        )

patient.rawdata.sg = rbind( patient.colon,
                            patient.rectum
                            )

dim( patient.rawdata.sg ) ## 253129     18
length( unique( patient.rawdata.sg$no ) ) ## 453 samples
length( unique( patient.rawdata.sg$item.c ) ) ## 331 items

rm( patient.colon.2011_2013,
    patient.colon.2014,
    patient.colon.2015,
    patient.colon.2016,
    patient.colon,
    patient.rectum.2011_2013,
    patient.rectum.2014,
    patient.rectum.2015,
    patient.rectum.2016,
    patient.rectum
    )


## ##############################################################################################################
################################ (2) extract the useful columns
head( patient.rawdata.sg )

patient.data.sg = patient.rawdata.sg[, c( "no",
                                          "location",
                                          "age",
                                          "sex", 
                                          "start_date_time",
                                          "stage_T","stage_N","stage_M","stages",
                                          "result_date_time", 
                                          "item.c",
                                          "value", 
                                          "ref"
                                          )
                                     ]

head(patient.data.sg)
dim( patient.data.sg ) ## 253129     13
length( unique( patient.data.sg$no ) ) ## 453 number of samples
length( unique( patient.data.sg$item.c ) )  ## 331 items

## ##############################################################################################################
##################################### (3) 校正分期
unique( patient.data.sg[, c("stage_T","stage_N","stage_M","stages")] )

patient.data.sg = patient.data.sg[-131878, ]

patient.data.sg$stages = "unknown"

#### 将NA保护起来
patient.data.sg$stage_M[is.na(patient.data.sg$stage_M)] = 1000
patient.data.sg$stage_N[is.na(patient.data.sg$stage_N)] = 1000
patient.data.sg$stage_T[is.na(patient.data.sg$stage_T)] = 1000

patient.data.sg[ patient.data.sg$stage_M == "1", "stages"] = "4"
patient.data.sg[ (patient.data.sg$stage_M == "0") &
                   ( (patient.data.sg$stage_N == "1")|
                     (patient.data.sg$stage_N == "2")|
                     (patient.data.sg$stage_N == "3")
                    ),
                 "stages"] = "3"

patient.data.sg[ (patient.data.sg$stage_M == "0") &
                 (patient.data.sg$stage_N == "0") &
                 ( (patient.data.sg$stage_T == "3") |
                   (patient.data.sg$stage_T == "4")
                   ),
                 "stages"] = "2"
patient.data.sg[ (patient.data.sg$stage_M == "0") &
                 (patient.data.sg$stage_N == "0") &
                 ( (patient.data.sg$stage_T == "1") |
                   (patient.data.sg$stage_T == "2")
                   ),
                 "stages"] = "1"
patient.data.sg[ (patient.data.sg$stage_M == "0") &
                  (patient.data.sg$stage_N == "0") &
                  (patient.data.sg$stage_T == "0"),
                 "stages"] = "NoTumor"

#### 将NA释放起来
patient.data.sg$stage_M[patient.data.sg$stage_M==1000] = NA
patient.data.sg$stage_N[patient.data.sg$stage_N==1000] = NA
patient.data.sg$stage_T[patient.data.sg$stage_T==1000] = NA

###  删除T,N,M 
patient.data.sg = patient.data.sg[, -c(6,7,8)]

head(patient.data.sg)
dim( patient.data.sg ) ## 253128     10
length( unique( patient.data.sg$no ) ) ## 452 number of samples
length( unique( patient.data.sg$item.c ) )  ## 330 items

## ############################################################################### 
##### (4) 筛选第一次手术之前的诊断化验

##### result_data_time出现了1900,视为错误删除
length( grep("1900",patient.data.sg$result_date_time) ) ## 3044
patient.data.sg$result_date_time[grep("1900",patient.data.sg$result_date_time)] = ""

###删除没有日期信息的
patient.data.sg = patient.data.sg[ (patient.data.sg$start_date_time != "") &
                             (patient.data.sg$result_date_time != "") , ]

dim( patient.data.sg ) ## 244646      10
length( unique( patient.data.sg$no) ) ## 442 number of samples
length( unique( patient.data.sg$item) )  ## 325 items

###筛选手术之前的检验
patient.data.sg$start_date_time  = as.POSIXct(patient.data.sg$start_date_time)

patient.data.sg = patient.data.sg[ order(patient.data.sg$result_date_time, decreasing = FALSE), ]

patient.data.sg01 = patient.data.sg[1:57321,]
patient.data.sg02 = patient.data.sg[57322:244646,]

patient.data.sg01$result_date_time = as.POSIXct(patient.data.sg01$result_date_time)
patient.data.sg02$result_date_time = as.POSIXct(patient.data.sg02$result_date_time)

patient.data.sg = rbind(patient.data.sg01, patient.data.sg02)

patient.data.sg = patient.data.sg[patient.data.sg$start_date_time >
                                    patient.data.sg$result_date_time, ]

dim( patient.data.sg ) ## 83442     10
length( unique( patient.data.sg$no) ) ## 436 number of samples
length( unique( patient.data.sg$item.c) )  ## 296 items


###筛选第一次手术之前的检验

###排序
patient.data.sg = patient.data.sg[ order( patient.data.sg$no,
                                    patient.data.sg$start_date_time,
                                    decreasing = FALSE
                                    ), 
                             ]

no              = unique(patient.data.sg$no)
id              = match( no, patient.data.sg$no )
start_date_time = patient.data.sg$start_date_time[id]

no.start.time.paste = data.frame( no = no, start.time = start_date_time )
tmp = Reduce( rbind,
              apply( no.start.time.paste, MARGIN = 1,
                     FUN = function(x){
                       return(
                         patient.data.sg[ patient.data.sg$no == x["no"] &
                                       patient.data.sg$start_date_time == x["start.time"],
                                       ]
                             )
                       }
                     )
              )

patient.data.sg = tmp[, c("no", "location", "stages","sex", "age", "item.c", "value", "ref" ) ]

head(patient.data.sg)
dim( patient.data.sg ) ## 77107     7
length( unique( patient.data.sg$no) ) ## 436 number of samples
length( unique( patient.data.sg$item) )  ## 296 items
rm(tmp)
rm(patient.data.sg01, patient.data.sg02)

## #####################################################################################
### (5) examine each columns
head( table(patient.data.sg$no) )
tail( table(patient.data.sg$no) )
sum( table(patient.data.sg$no) )
length( unique( patient.data.sg$no) ) ## 436 number of samples

table(patient.data.sg$location)
sum(table(patient.data.sg$location))
sum(patient.data.sg$location == "") ## 0

table( patient.data.sg$sex )
sum(table( patient.data.sg$sex ))

table( patient.data.sg$age )
sum(table( patient.data.sg$age ))

########################################## item
unique( patient.data.sg$item.c ) ## 296 items
length(unique( patient.data.sg$item.c ))
sum(table( patient.data.sg$item.c ))
sum( patient.data.sg$item.c == "") ## 10
patient.data.sg$item.c[ patient.data.sg$item.c == "" ] = NA

# ## delete "*" or " " of start of some items
unique( patient.data.sg$item.c[ grep( "^[*]", patient.data.sg$item.c ) ] )
## 3个指标分带和不带“*”：红细胞，白细胞，葡萄糖
table(patient.rawdata.sg[patient.rawdata.sg$item.c == "*红细胞","specimen"]) ## 血
table(patient.rawdata.sg[patient.rawdata.sg$item.c == "红细胞","specimen"])  ## 尿
table(patient.rawdata.sg[patient.rawdata.sg$item.c == "*白细胞","specimen"]) ## 血
table(patient.rawdata.sg[patient.rawdata.sg$item.c == "白细胞","specimen"])  ## 尿
table(patient.rawdata.sg[patient.rawdata.sg$item.c == "* 葡萄糖","specimen"]) ## 血清
table(patient.rawdata.sg[patient.rawdata.sg$item.c == "葡萄糖","specimen"])   ## 全血

patient.data.sg$item.c[patient.data.sg$item.c == "红细胞"]   = "尿红细胞"
patient.data.sg$item.c[patient.data.sg$item.c == "白细胞"]   = "尿白细胞"
patient.data.sg$item.c[patient.data.sg$item.c == "* 葡萄糖"] = "葡萄糖"

patient.data.sg$item.c[ grep( "^[*]", patient.data.sg$item.c ) ] =
  gsub( "^[*]", "", patient.data.sg$item.c[grep("^[*]", patient.data.sg$item.c)] )

patient.data.sg$item.c[ grep( "^[ ]", patient.data.sg$item.c ) ] =
  gsub( "^[ ]", "", patient.data.sg$item.c[grep("^[ ]", patient.data.sg$item.c)] )


unique( patient.data.sg$item.c ) ## 295 items
sum(table( patient.data.sg$item.c ))

#################### value
head(table( patient.data.sg$value ),10)
tail(table( patient.data.sg$value ),10)
sum( table( patient.data.sg$value ) ) ## 77107
sum( is.na( patient.data.sg$value ) ) ## 0

patient.data.sg$value[patient.data.sg$value == 1000]
patient.data.sg$value[patient.data.sg$value == 2000]

##  先把阴阳保护起来，可以继续推敲
patient.data.sg$value[ grep( "[阴]", patient.data.sg$value ) ] = 1000
patient.data.sg$value[ grep( "[阳]", patient.data.sg$value ) ] = 2000

patient.data.sg$value[ patient.data.sg$value == "-----" ]     = 1000
patient.data.sg$value[ patient.data.sg$value == "---" ]       = 1000
patient.data.sg$value[ patient.data.sg$value == "-性" ]       = 1000
patient.data.sg$value[ grep( "[+]", patient.data.sg$value ) ] = 2000

patient.data.sg$value[ patient.data.sg$value == "NEGATIVE" ]  = 1000
patient.data.sg$value[ patient.data.sg$value == "POSITIVE" ]  = 1000

patient.data.sg$value = as.numeric( patient.data.sg$value )

patient.data.sg$value[ patient.data.sg$value == 1000 ] = "阴"
patient.data.sg$value[ patient.data.sg$value == 2000 ] = "阳"
## 问题：数值型指标的value有的是用阴阳标记的，这部分也被保护了
## 解决：最后将将数值型指标再as.numeric(),阴阳变成了NA 

#################### ref
table( patient.data.sg$ref )
patient.data.sg$ref[ grep( "[阴]", patient.data.sg$ref ) ] = "阴性"

### delete the data,whose "no","sex","age",or "value" are empty
patient.data.sg = na.omit( patient.data.sg, cols = c("no", "sex", "age", "item.c", "value") )

head(patient.data.sg)
dim( patient.data.sg ) ## 73268     8
length( unique( patient.data.sg$no) ) ## 436 number of samples
length( unique( patient.data.sg$item.c) )  ## 280 items

####### delete repeat data,just reserve the newest one
head(patient.data.sg)
patient.data.sg.paste = paste0( patient.data.sg$no, "_", 
                                patient.data.sg$item
                                )
id                 = match( unique(patient.data.sg.paste), patient.data.sg.paste )
patient.data.sg       = patient.data.sg[ id, ]

dim( patient.data.sg ) ## 44966     8
length( unique( patient.data.sg$no) ) ## 436 number of samples
length( unique( patient.data.sg$item.c) )  ## 280 items
range( table( patient.data.sg$no ) )  ## 18 171 range of items

rownames(patient.data.sg) = NULL

write.table( patient.data.sg,
             "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/patient.data.sg.txt",
             quote = T )

## #########################################################################
####### (6) 结合patient与health数据确定研究指标
length( unique( patient.data.sg$no) ) ## 436 number of samples
range( table( patient.data.sg$no ) )  ## 18 171 range of items
length( table(patient.data.sg$item.c) ) ## 280 items

sort( table( patient.data.sg$item.c ), decreasing=TRUE )
# sum( table(patient.data.sg$item.c) >= 340 ) ## 71 items whose samples are more than 100

patient.item.35.sg = rownames( 
  sort(
    table(patient.data.sg$item),decreasing=TRUE
  )
)[1: floor((length(unique(patient.data.sg$item)))*0.35)]

write.table( patient.item.35.sg,
             "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/patient.item.35.sg.txt",
             quote = T )

health.item.40.sg = read.table(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/health/health.item.40.sg.txt",
  stringsAsFactors = F
  ) ## 80 items

health.item.40.sg = c(health.item.40.sg$x)

items.sg = intersect( patient.item.35.sg, health.item.40.sg ) ## 40 items

write.table( items.sg,
             "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/items.sg.txt",
             quote = T )

######################################################################################
## (7) 筛选 28 items 的数据 
patient.data.sg.new = as.data.frame(matrix(numeric(0),ncol=4)) ## data about 41 items
patient.ref.sg = as.data.frame(matrix(numeric(0),ncol=2))
for(i in items.sg){
  patient.data.sg.new.single = patient.data.sg[ patient.data.sg$item == i, ]
  patient.ref.sg             = rbind(patient.ref.sg, patient.data.sg.new.single[1,c("item.c","ref")])
  patient.data.sg.new        = rbind( patient.data.sg.new, patient.data.sg.new.single )
}

patient.data.sg.new = patient.data.sg.new[order(patient.data.sg.new$no), ]
rownames(patient.data.sg.new) = NULL
rownames(patient.ref.sg)      = NULL
rm( patient.data.sg.new.single )

dim(patient.data.sg.new) ## 11658     8
length( unique( patient.data.sg.new$no ) ) ## 436 number of samples
range( table( patient.data.sg.new$no ) )  ## 2 28 range of items
length( unique( patient.data.sg.new$item.c) )  ## 28 items

write.table(patient.ref.sg, 
            "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/patient.ref.sg.txt",
            quote = T, sep = "\t"
            )

write.table(patient.data.sg.new, 
            "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/patient.data.sg.new.txt",
            quote = F
            )

############################################################################
## (8) transform raw data "patient.data.sg_new" into "patient.matrix"
library(reshape2)
patient.info.sg = dcast(patient.data.sg.new, no+location+stages+sex+age  ~ item.c, value = "value")

###  查看指标是定性还是定量

#### 将变量值由character变成numeric
items.num = setdiff( colnames(patient.info.sg), c("no", "location", "stages", "sex") )
for (i in items.num) {
  patient.info.sg[, i] = as.numeric( patient.info.sg[, i] )
  }

write.table(patient.info.sg, 
            "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/patient.info.sg.txt", 
            quote = T
            )

############################################################################
save.image(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/sg_hospital_patient(10-15).RData"
  )


