# preprocess of patient data of zl hospital 2010-2015
## 2016.04, shiq

setwd("/lustre/user/liclab/shiq/CRC/zl_hospital/patient")

## #################################################################################
####### (1) import data 2010-2015
patient.2010 = read.csv(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/CRC_2010-2015/2010出院的肠含诊断-his.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
  )
patient.2011 = read.csv(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/CRC_2010-2015/2011出院的肠含诊断-his.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.2012 = read.csv(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/CRC_2010-2015/2012出院的肠含诊断-his.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.2013 = read.csv(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/CRC_2010-2015/2013出院的结肠含诊断ID加密非术后-his.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.2014 = read.csv(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/CRC_2010-2015/2014出院的结肠含诊断ID加密非术后-his.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
patient.2015 = read.csv(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/CRC_2010-2015/2015出院的结肠含诊断ID加密非术后-his.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)

patient.rawdata.zl = rbind( patient.2010,
                         patient.2011,
                         patient.2012,
                         patient.2013,
                         patient.2014,
                         patient.2015
                         )
dim( patient.rawdata.zl ) ## 1354200      18
length( unique( patient.rawdata.zl$patient_id ) ) ## 7068 samples
length( unique( patient.rawdata.zl$report_item_name ) ) ## 363 items

rm( patient.2010,
    patient.2011,
    patient.2012,
    patient.2013,
    patient.2014,
    patient.2015
    )

## ##############################################################################################################
####### (2) screening stage data with tumor_stage_other(综合分期)

## 分期信息 importing
patient.stage.raw = read.csv(
  "/lustre/user/liclab/shiq/CRC/sg_hospital/patient/2010-2015出院含肿瘤分期并脱敏.txt",
  header = T,sep = "\t", fill = T, stringsAsFactors = F
)
length( unique( patient.stage.raw$patient_id ) ) ## 5329 samples

## 分期信息with tumor_stage_other(综合分期)
sum(table(patient.stage.raw$patient_id))
patient.stage = patient.stage.raw[, c( "patient_id", "visit_id", "tumor_stage_other" )]

## 删除重复的的分期信息, 完整的分期信息
patient.stage = patient.stage[ order(patient.stage$visit_id,decreasing=FALSE), ]
patient.stage = patient.stage[ !duplicated( patient.stage$patient_id ), ]
dim( patient.stage ) ## 5329    3

table(patient.stage$tumor_stage_other)
patient.stage$tumor_stage_other[ patient.stage$tumor_stage_other == "Ⅲ期 "]  =  "Ⅲ期"
patient.stage$tumor_stage_other[ patient.stage$tumor_stage_other == "Ⅱ期 "]  =  "Ⅱ期"


## ##############################################################################################################
##### (3) extract the useful columns
head( patient.rawdata.zl )

patient.data.zl = patient.rawdata.zl[, c( "patient_id", 
                                       "diagnosis_desc", 
                                       "age",
                                       "sex", 
                                       "start_date_time",
                                       "result_date_time",
                                       "report_item_name",
                                       "result",
                                       "normal_value_interval" 
                                       )
                                  ]
colnames(patient.data.zl) = c( "no", "location","age", "sex", "start_date_time",
                            "result_date_time", "item", "value", "ref"
                            )
head(patient.data.zl)
dim( patient.data.zl ) ## 1354200       9
length( unique( patient.data.zl$no ) ) ## 7068 number of samples
length( unique( patient.data.zl$item ) )  ## 363 items

## ############################################################################### 
##### (4) 筛选第一次手术之前的诊断化验

###删除没有日期信息的
patient.data.zl = patient.data.zl[ patient.data.zl$start_date_time != "" &
                             patient.data.zl$result_date_time != "" , ]

dim( patient.data.zl ) ## 1156152       9
length( unique( patient.data.zl$no) ) ## 5339 number of samples
length( unique( patient.data.zl$item) )  ## 337 items

###筛选手术之前的检验
patient.data.zl$start_date_time  = as.POSIXct(patient.data.zl$start_date_time)
patient.data.zl$result_date_time = as.POSIXct(patient.data.zl$result_date_time)

patient.data.zl = patient.data.zl[patient.data.zl$start_date_time >
                            patient.data.zl$result_date_time, ]

dim( patient.data.zl ) ## 406599      9
length( unique( patient.data.zl$no) ) ## 4380 number of samples
length( unique( patient.data.zl$item) )  ## 305 items

###筛选第一次手术之前的检验

###排序
patient.data.zl = patient.data.zl[ order( patient.data.zl$no,
                                    patient.data.zl$start_date_time,
                                    decreasing = FALSE
                                    ), 
                             ]

no              = unique(patient.data.zl$no)
id              = match( no, patient.data.zl$no )
start_date_time = patient.data.zl$start_date_time[id]

no.start.time.paste = data.frame( no = no, start.time = start_date_time )
tmp = Reduce( rbind,
              apply( no.start.time.paste, MARGIN = 1,
                     FUN = function(x){
                       return(
                         patient.data.zl[ patient.data.zl$no == x["no"] &
                                       patient.data.zl$start_date_time == x["start.time"],
                                       ]
                             )
                       }
                     )
              )

patient.data.zl = tmp[, c("no", "location", "sex", "age", "item", "value", "ref" ) ]

head(patient.data.zl)
dim( patient.data.zl ) ## 358797      7
length( unique( patient.data.zl$no) ) ## 4380 number of samples
length( unique( patient.data.zl$item) )  ## 292 items
rm(tmp)

## #####################################################################################
### (5) examine each columns
head( table(patient.data.zl$no) )
tail( table(patient.data.zl$no) )
sum( table(patient.data.zl$no) )
length( unique( patient.data.zl$no) ) ## 4380 number of samples

sum(table(patient.data.zl$location))
sum(patient.data.zl$location == "") ## 0

table( patient.data.zl$sex )
sum(table( patient.data.zl$sex ))

table( patient.data.zl$age )
sum(table( patient.data.zl$age ))

#################### item
unique( patient.data.zl$item ) ## 292 items
length(unique( patient.data.zl$item ))
sum(table( patient.data.zl$item ))
sum( patient.data.zl$item == "") ## 36
patient.data.zl$item[ patient.data.zl$item == "" ] = NA

## delete "." of tail of some items
patient.data.zl$item[ grep( "[.]$", patient.data.zl$item ) ] = 
  gsub( "[.]$", "", patient.data.zl$item[grep("[.]$", patient.data.zl$item)] )
patient.data.zl$item[ patient.data.zl$item == "T-Bil" ] = "TBil"
patient.data.zl$item[ patient.data.zl$item == "PH" ] = "pH"

unique( patient.data.zl$item ) ## 260 items
sum(table( patient.data.zl$item ))

#################### value
head(table( patient.data.zl$value ))
tail(table( patient.data.zl$value ))
sum( table( patient.data.zl$value ) ) ## 358787
sum( is.na( patient.data.zl$value ) ) ## 0

patient.data.zl$value[patient.data.zl$value == 1000]
patient.data.zl$value[patient.data.zl$value == 2000]

##  先把阴阳保护起来，可以继续推敲
patient.data.zl$value[ grep( "[阴]", patient.data.zl$value ) ] = 1000
patient.data.zl$value[ grep( "[阳]", patient.data.zl$value ) ] = 2000

patient.data.zl$value[ patient.data.zl$value == "-----" ]     = 1000
patient.data.zl$value[ patient.data.zl$value == "---" ]       = 1000
patient.data.zl$value[ patient.data.zl$value == "-性" ]       = 1000
patient.data.zl$value[ grep( "[+]", patient.data.zl$value ) ] = 2000

patient.data.zl$value[ patient.data.zl$value == "NEGATIVE" ]  = 1000
patient.data.zl$value[ patient.data.zl$value == "POSITIVE" ]  = 1000

patient.data.zl$value = as.numeric( patient.data.zl$value )

patient.data.zl$value[ patient.data.zl$value == 1000 ] = "阴"
patient.data.zl$value[ patient.data.zl$value == 2000 ] = "阳"
## 问题：数值型指标的value有的是用阴阳标记的，这部分也被保护了
## 解决：最后将将数值型指标再as,numeric(),阴阳变成了NA 

#################### ref
table( patient.data.zl$ref )
patient.data.zl$ref[ grep( "[阴]", patient.data.zl$ref ) ] = "阴性"

### delete the data,whose "no","sex","age",or "value" are empty
patient.data.zl = na.omit( patient.data.zl, cols = c("no", "sex", "age", "item", "value") )

head(patient.data.zl)
dim( patient.data.zl ) ## 353304      7
length( unique( patient.data.zl$no) ) ## 4378 number of samples
length( unique( patient.data.zl$item) )  ## 249 items

####### delete repeat data,just reserve the newest one
head(patient.data.zl)
patient.data.zl.paste = paste0( patient.data.zl$no, "_", 
                                patient.data.zl$item
                                )
id                 = match( unique(patient.data.zl.paste), patient.data.zl.paste )
patient.data.zl       = patient.data.zl[ id, ]

dim( patient.data.zl ) ## 307778      7
length( unique( patient.data.zl$no) ) ## 4378 number of samples
length( unique( patient.data.zl$item) )  ## 249 items
range( table( patient.data.zl$no ) )  ## 1 117 range of items

rownames(patient.data.zl) = NULL

write.table( patient.data.zl,
             "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/patient.data.zl.txt",
             quote = T )

## #########################################################################
####### (6) 结合patient与health数据确定研究指标
length( unique( patient.data.zl$no) ) ## 4378 number of samples
range( table( patient.data.zl$no ) )  ## 1 119 range of items
length( unique(patient.data.zl$item) ) ## 251 items

sort( table( patient.data.zl$item ), decreasing=TRUE )
# sum( table(patient.data.zl$item) >= 950 ) ## 88 items whose samples are more than 950

patient.item.35.zl = rownames( 
  sort(
    table(patient.data.zl$item),decreasing=TRUE
    )
  )[1: floor((length(unique(patient.data.zl$item)))*0.35)]

write.table( patient.item.35.zl,
             "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/patient.item.35.zl.txt",
             quote = T )

health.item.40.zl = read.table(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/health/health.item.40.zl.txt",
  stringsAsFactors = F
  ) ## 67 items

health.item.40.zl = c(health.item.40.zl$x)

items.zl = intersect( patient.item.35.zl, health.item.40.zl ) ## 40 items

write.table( items.zl,
             "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/items.zl.txt",
             quote = T )


######################################################################################
## (7) 筛选 40 items 的数据 
patient.data.zl.new = as.data.frame(matrix(numeric(0),ncol=4)) ## data about 41 items
patient.ref.zl = as.data.frame(matrix(numeric(0),ncol=2))
for(i in 1:length(items.zl)){
  patient.data.zl.new.single = patient.data.zl[ patient.data.zl$item == items.zl[i], ]
  patient.ref.zl             = rbind(patient.ref.zl, patient.data.zl.new.single[1,c("item","ref")])
  patient.data.zl.new        = rbind( patient.data.zl.new, patient.data.zl.new.single )
}

patient.data.zl.new = patient.data.zl.new[order(patient.data.zl.new$no), ]
rownames(patient.data.zl.new) = NULL
rownames(patient.ref.zl)      = NULL
rm( patient.data.zl.new.single )

dim(patient.data.zl.new) ## 144311      7
length( unique( patient.data.zl.new$no ) ) ## 4211 number of samples
range( table( patient.data.zl.new$no ) )  ## 1 40 range of items
length( unique( patient.data.zl.new$item) )  ## 40 items

write.table(patient.ref.zl, 
            "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/patient.ref.zl.txt",
            quote = T, sep = "\t"
            )

write.table(patient.data.zl.new, 
            "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/patient.data.zl.new.txt",
            quote = F
            )

############################################################################
## (8) transform raw data "patient.data.zl_new" into "patient.matrix"
library(reshape2)
tmp = dcast(patient.data.zl.new, no + location  ~ item, value = "value")

####### 添加分期信息,correct the wrong age and sex
id    = match( tmp$no, patient.stage$patient_id )
stage = patient.stage$tumor_stage_other[id]

id    = match( tmp$no, patient.data.zl.new$no )
sex   = patient.data.zl.new$sex[id]
age   = patient.data.zl.new$age[id]

head(tmp)
patient.info.zl = cbind( tmp[, c("no", "location")],
                         stage,
                         sex,
                         age,
                         tmp[, -c(1,2)], 
                         stringsAsFactors = FALSE
                         )
head(patient.info.zl)
patient.info.zl$stage[ is.na(patient.info.zl$stage) ] = "unknown"

###  查看指标是定性还是定量
###  HBeAg HBsAg在patient中是定量，而在health中基本全是定性，且ref不一致，故删除
patient.info.zl = patient.info.zl[, setdiff( colnames(patient.info.zl), c("HBeAg", "HBsAg")) ]

#### 将变量值由character变成numeric
items.num = setdiff( colnames(patient.info.zl), c("no", "location", "stage", "sex") )
for (i in items.num) {
  patient.info.zl[, i] = as.numeric( patient.info.zl[, i] )
}


write.table(patient.info.zl, 
            "/lustre/user/liclab/shiq/CRC/zl_hospital/patient/patient.info.zl.txt", 
            quote = T
            )
rm( tmp )

############################################################################
save.image(
  "/lustre/user/liclab/shiq/CRC/zl_hospital/zl_hospital_patient(10-15).RData"
  )


