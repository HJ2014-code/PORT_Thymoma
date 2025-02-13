rm(list=ls())
# 设置环境参数
work_dir <- "C:/Users/huangjing/Desktop/01_规培/科研/CHART_PORT"
# 设置工作目录
setwd(work_dir)


thymus <- read.csv("01_Data.csv",check.names=FALSE)

table(thymus$`最终切除状况(0-R0,1-1,2-2,3-biopsy)`)
thymoma_R0 <- thymus[thymus$`最终切除状况(0-R0,1-1,2-2,3-biopsy)` == "0",]
write.csv(thymoma_R0, file = "01_Data_R0.csv")

#####
#载入所需的R包
library("survival")
library("survminer")
#载入并查看数据集
thymus_ <- read.csv("03_Data_full_info_1853.csv",check.names=FALSE)
thymus <- thymus_[thymus_$OS_month != "0",]

table(thymus$Masaoka)
thymoma_II <- thymus[thymus$Masaoka == "II",]

attach(thymoma_II)
fit <- survfit(Surv(OS_month,OS_STATUS) ~ Radio,
               data=thymoma_II)

plot(fit)
ggsurvplot(fit)
#优化
p3<-ggsurvplot(fit, data = thymoma_II,
               pval = TRUE,
               pval.coord = c(0.5,0.05),
               risk.table = TRUE,
               xlab = "Follow up time (months)",
               break.x.by = 24,
               legend = c(0.8,0.75),
               legend.title = "",
               legend.labs = c("Surgery without PORT","Surgery with PORT"),
               palette="lancet", ylab = "Overall survival" )
p3
###添加COX回归hazard ratio值相关信息
res_cox<-coxph(Surv(OS_month,OS_STATUS) ~ Radio, 
               data=thymoma_II)
p3$plot = p3$plot + ggplot2::annotate("text",x = 50, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 50, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 50, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3


rm(list=ls())
thymoma_II <- read.csv("05_thymoma_II.csv",check.names=FALSE)
library(MatchIt)
###开始进行匹配
###先将data中的group定义为逻辑型变量，
#如果为”uncensored”则返回’TRUE’，否则返回’FALSE’

###这里一定注意定义为’TRUE’的变量必须少于定义为’FALSE’的变量
### y~x分别为干预（分组）因素和需要匹配的变量，
#即group2为干预（分组）因素，gender、age和TNM则作为需要匹配的变量，
#然后将文件赋值为data_match

### method为匹配算法，具体参见MatchIt包的说明，一般采用“近邻法”进行匹配

### ratio为匹配的比例，默认为1:1，一般来说也这样选择。
#如果对照组人数是干预组的若干倍，那ratio的值可以稍大一点，如2:1
###caliper代表的是卡钳值，反应实验组与对照组在进行配对时允许的误差，
#如果设定为0.02，则实验组与对照组按照倾向性评分±0.02进行匹配
#（如果不加卡钳值，则实验组全部被匹配）
table(thymoma_II$Radio)
thymoma_II$Radio2 <- as.logical(thymoma_II$Radio == '1')
table(thymoma_II$Radio2)
write.csv(thymoma_II, file = "05_thymoma_II.csv")

median(thymoma_II$Age_primary)
table(thymoma_II$WHO_classification_Combine)

attach(thymoma_II)
data_match <- matchit(Radio2~Age+Sex+Surg_Extent+WHO_classification_Combine
                      +Tumor_Size+Chemotherapy, 
                      data = thymoma_II, 
                      method='nearest', ratio=2, caliper = 0.02)
summary(data_match) #汇总结果，重点关注红框里的内容


dta_m <- match.data(data_match)
###KM thymus
fit <- survfit(Surv(OS_month,OS_STATUS) ~ Radio, 
               data=dta_m)
plot(fit)
ggsurvplot(fit)

#优化
p3<-ggsurvplot(fit, data = dta_m,
               pval = TRUE,
               pval.coord = c(0.5,0.05),
               risk.table = TRUE,
               xlab = "Follow up time (months)",
               break.x.by = 24,
               legend = c(0.8,0.75),
               legend.title = "",
               legend.labs = c("Surgery without PORT","Surgery with PORT"),
               palette="lancet", ylab = "Overall survival" )
p3
###添加COX回归hazard ratio值相关信息
res_cox<-coxph(Surv(OS_month,OS_STATUS) ~ Radio, 
               data=dta_m)
p3$plot = p3$plot + ggplot2::annotate("text",x = 50, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 50, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 50, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3

####Baseline Characteristics#######
###################################
###################################
library(tableone)
## 需要汇总的变量Age+Sex+Surg_Extent+WHO_classification_Combine
#+Tumor_Size+Chemotherapy
myVars <- c('Age','Sex', 'Surg_Extent','WHO_classification_Combine',
            'Tumor_Size','Chemotherapy')
## 需要转为分类变量的变量
catVars <- c('Age','Sex', 'Surg_Extent','WHO_classification_Combine',
             'Tumor_Size','Chemotherapy')
## Create a TableOne object
tab2 <- CreateTableOne(vars = myVars, 
                       data = dta_m, 
                       factorVars = catVars,
                       strata = 'Radio',
                       addOverall = TRUE)

tab2 <- CreateTableOne(vars = myVars, 
                       data = thymoma_II, 
                       factorVars = catVars,
                       strata = 'Radio',
                       addOverall = TRUE)
tab2
###数据导出
tab2Mat <- print(tab2,quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE)
## 保存为 CSV 格式文件
write.csv(tab2Mat, file = "03_myTable_thymoma_Baseline_PSM.csv")
write.csv(tab2Mat, file = "03_myTable_thymoma_Baseline_before_PSM.csv")



#####单因素分析####
###################
###################
rm(list=ls())
thymoma_II <- read.csv("05_thymoma_II.csv",check.names=FALSE)

#1-1.批量单因素Cox包
library(survival)
library(plyr)

#单因素
# Age+Sex+Surg_Extent+WHO_classification_Combine
#+Tumor_Size+Chemotherapy
attach(thymoma_II)
table(thymoma_II$Tumor_Size_cm)
res.cox1 <- coxph(Surv(OS_month,OS_STATUS) ~ Tumor_Size, data = thymoma_II)
res.cox1
summary(res.cox1)

covariates <-c("Age","Sex","Surg_Extent", 
               "WHO_classification_Combine","Tumor_Size", "Chemotherapy","Radio")

## 得到一个列表分别为，对每个变量构建的生存对象公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS_month,OS_STATUS)~', x)))
#univ_formulas
## 巧妙的配合lapply                        
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = thymoma_II)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-round(x$wald["pvalue"], digits=4)
                         wald.test<-round(x$wald["test"], digits=4)
                         beta<-round(x$coef[1], digits=4);#coeficient beta
                         HR <-round(x$coef[2], digits=4);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],4)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
                         HR_ <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR_, HR, HR.confint.lower,HR.confint.upper,wald.test, p.value)
                         names(res)<-c("beta", " HR_","HR","HR.confint.lower","HR.confint.upper","wald.test", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

class(univ_results)
str(univ_results)
res <- t(as.data.frame(univ_results, check.names = FALSE))
res<-as.data.frame(res)
res <- res[order(res$p.value),]
write.csv(res,file = "05_uni_results.csv")


#####subgroup分析####
###################
###################
rm(list=ls())
thymoma_II <- read.csv("05_thymoma_II.csv",check.names=FALSE)

#变量Age+Sex+Surg_Extent+WHO_classification_Combine
#+Tumor_Size+Chemotherapy+stage IIA  IIB
table(thymoma_II$Masaoka_Subgroup)
table(thymoma_II$Age)
table(thymoma_II$Sex)
table(thymoma_II$Surg_Extent)
table(thymoma_II$WHO_classification_Combine)
table(thymoma_II$WHO_classification_Combine_B)
table(thymoma_II$Tumor_Size)
table(thymoma_II$Chemotherapy)


IIA <- thymoma_II[thymoma_II$Masaoka_Subgroup == "IIA",]
IIA <- thymoma_II[thymoma_II$Masaoka_Subgroup == "IIB",]
IIA <- thymoma_II[thymoma_II$Masaoka_Subgroup == "II",]

IIA <- thymoma_II[thymoma_II$Age == "<60",]
IIA <- thymoma_II[thymoma_II$Age == "≥60",]

IIA <- thymoma_II[thymoma_II$Sex == "1",]
IIA <- thymoma_II[thymoma_II$Sex == "2",]

IIA <- thymoma_II[thymoma_II$Surg_Extent == "1",]
IIA <- thymoma_II[thymoma_II$Surg_Extent == "2",]

IIA <- thymoma_II[thymoma_II$WHO_classification_Combine == "A+AB",]
IIA <- thymoma_II[thymoma_II$WHO_classification_Combine == "B1+B2+B3",]

IIA <- thymoma_II[thymoma_II$WHO_classification_Combine_B == "A+AB+B1",]
IIA <- thymoma_II[thymoma_II$WHO_classification_Combine_B == "B2+B3",]

IIA <- thymoma_II[thymoma_II$Tumor_Size == "<6.5",]
IIA <- thymoma_II[thymoma_II$Tumor_Size == "≥6.5",]

IIA <- thymoma_II[thymoma_II$Chemotherapy == "0",]
IIA <- thymoma_II[thymoma_II$Chemotherapy == "2",]

#载入所需的R包
library("survival")
library("survminer")
attach(IIA)
fit <- survfit(Surv(OS_month,OS_STATUS) ~ Radio,
               data=IIA)

plot(fit)
ggsurvplot(fit)
#优化
p3<-ggsurvplot(fit, data = IIA,
               pval = TRUE,
               pval.coord = c(0.5,0.05),
               risk.table = TRUE,
               xlab = "Follow up time (months)",
               break.x.by = 24,
               legend = c(0.8,0.75),
               legend.title = "",
               legend.labs = c("Surgery without PORT","Surgery with PORT"),
               palette="lancet", ylab = "Overall survival" )
p3
###添加COX回归hazard ratio值相关信息
res_cox<-coxph(Surv(OS_month,OS_STATUS) ~ Radio, 
               data=IIA)
p3$plot = p3$plot + ggplot2::annotate("text",x = 50, y = 0.15,
                                      label = paste("HR :",round(summary(res_cox)$conf.int[1],2))) + ggplot2::annotate("text",x = 50, y = 0.10,
                                                                                                                       label = paste("(","95%CI:",round(summary(res_cox)$conf.int[3],2),"-",round(summary(res_cox)$conf.int[4],2),")",sep = ""))+
  ggplot2::annotate("text",x = 50, y = 0.05,
                    label = paste("P:",round(summary(res_cox)$coef[5],4)))
p3



####subgroup 结果可视化
rm(list=ls())
dataxx <-read.csv("06_Subtype_P.csv",check.names=FALSE)
datax <- dataxx [,c(4)]
datax <- as.matrix(dataxx [,4])
barplot(datax,beside = TRUE,horiz =T)
barplot(datax,beside = TRUE,angle = 15+10*1:5, density = 10,col = rainbow(2),horiz =T)
barplot(datax,beside = TRUE,
        xlim = c(0,1.5), 
        #names.arg = dataxx [,2],  # 柱子名称
        angle = 15+10*1:5, density = 10,,horiz =T)










