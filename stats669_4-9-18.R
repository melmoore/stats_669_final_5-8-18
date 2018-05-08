#Pres for stats 669 class 4-10-18



#load libraries

library(scales)
library(readr)
library(plyr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(tidyr)
library(reshape2)
library(cowplot)
library(nlme)
library(lme4)

library(extrafont)


#------------


#load data

tv <- read_csv("~/25-28-30_tv-final_clean.csv", 
               col_types = cols(temp.avg = col_factor(levels = c("25","28", "30")), 
                                temp.var = col_factor(levels = c("0", "5", "10")), 
                                treatment = col_factor(levels = c("control","para"))))

View(tv)



tv.long <- read_csv("~/25-28-30_tv-final_clean_LONG.csv", 
                    col_types = cols(temp.avg = col_factor(levels = c("25", "28", "30")), 
                                     temp.var = col_factor(levels = c("0","5", "10")), 
                                     treatment = col_factor(levels = c("control","para"))))

View(tv.long)



#Logging mass for plotting:

tv.long$log.mass<-log(tv.long$mass)



#---------



#Creating a few columns with treatment combinations to be able to sort out the 30+/-5 treatments

#Long dataset
tv.long<-unite(tv.long,tmp.trt,temp.avg,temp.var,sep=".",remove = FALSE)
tv.long<-unite(tv.long,all.trt,tmp.trt,treatment,sep=".",remove=FALSE)
tv.long$all.trt<-gsub(".para",".p",tv.long$all.trt)
tv.long$all.trt<-gsub(".control",".c",tv.long$all.trt)


tv.long.ab<-subset(tv.long,all.trt!="30.5.p")
tv.long.ab<-subset(tv.long.ab,all.trt!="30.5.c")


#Wide dataset
tv<-unite(tv,tmp.trt,temp.avg,temp.var,sep=".",remove = FALSE)
tv<-unite(tv,all.trt,tmp.trt,treatment,sep=".",remove=FALSE)
tv$all.trt<-gsub(".para",".p",tv$all.trt)
tv$all.trt<-gsub(".control",".c",tv$all.trt)


tv.ab<-subset(tv,all.trt!="30.5.p")
tv.ab<-subset(tv.ab,all.trt!="30.5.c")




#-------------

#plotting log.mass by age--includes cats without wasp emergence up to and including 5th instar--emergence points
  ##only include hosts with wasp emergence--no mongos

tv.long.ab<-subset(tv.long.ab,bug.id!="30.10_p_38")
tv.long.ab<-subset(tv.long.ab,bug.id!="30.10_p_23")
tv.long.ab<-subset(tv.long.ab,bug.id!="30.10_p_61")


amass.sum2<-summarySE(tv.long.ab,measurevar = "log.mass",
                      groupvar=c("treatment","temp.avg","temp.var","instar"),
                      na.rm=TRUE)
amass.sum2


#removing end day.age for those without wasp emergence at 30.10.p
tv.long.ab$suc.ovp[is.na(tv.long.ab$suc.ovp)]<-0
tv.long.ab$end.use<-ifelse(tv.long.ab$instar=="end" & tv.long.ab$all.trt=="30.10.p" & tv.long.ab$suc.ovp=="0",0,1)
tv.long.ab$suc.ovp[(tv.long.ab$suc.ovp)==0]<-NA
tv.long.ab<-subset(tv.long.ab,end.use==1)

dagem.sum2<-summarySE(tv.long.ab,measurevar = "day.age",
                      groupvar=c("treatment","temp.avg","temp.var","instar"),
                      na.rm=TRUE)
dagem.sum2


amass.sum2$day.age<-dagem.sum2[,6]
amass.sum2$dage.se<-dagem.sum2[,8]



amass.plot2<-ggplot(amass.sum2,aes(x=day.age,y=log.mass,
                                   group=interaction(temp.var,treatment),color=temp.var))
amass.plot2+geom_point(aes(shape=treatment),size=4
)+geom_line(aes(linetype=treatment),size=1.4
)+geom_errorbar(aes(ymin=log.mass-se,ymax=log.mass+se),width=.3,size=1
)+geom_errorbarh(aes(xmin=day.age-dage.se,xmax=day.age+dage.se)      
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10")
)+scale_linetype_manual(values=c("solid","dashed"),name="Treatment",
                        breaks=c("control","para"),labels=c("Cotrol","Parasitized")
)+scale_shape_manual(values = c(16,17),name="Treatment",
                     breaks=c("control","para"),labels=c("Cotrol","Parasitized")
)+labs(x="Age [days]",y="Log(Mass) [mg]"
)+facet_wrap(~temp.avg
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=24),
        axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.key.width=unit(15,"mm"),
        legend.background = element_rect(color="black",linetype="solid"))



#----------

#Plotting a frequency distribution curve of 30.10 dev time to wandering and emergence

all.30<-subset(tv,temp.avg=="30")
all.30<-subset(all.30,temp.var!="5")
all.30.c<-subset(all.30,treatment=="control")
all.30.p<-subset(all.30,treatment=="para")

#Removing those that wandered
all.30.p$mass.wand[is.na(all.30.p$mass.wand)]<-0
all.30.p<-subset(all.30.p,mass.wand==0)
all.30.p$mass.wand[all.30.p$mass.wand==0]<-NA

all.30.p$suc.ovp[is.na(all.30.p$suc.ovp)]<-1.5
all.30.p<-subset(all.30.p, suc.ovp>0)
all.30.p$suc.ovp[all.30.p$suc.ovp==1.5]<-NA


#Adding columns indicating whether wasps emerged ("y") or not ("n")
all.30.p$em.0<-ifelse(all.30.p$num.em=="0","n","y")
all.30.p$load.0<-ifelse(all.30.p$load=="0","n","y")

all.30.p$em.0<-factor(all.30.p$em.0,levels=c("y","n"))



#Fd plot colored by combination of temp.var and load.0

all.30.p<-unite(all.30.p,trt.load,temp.var,load.0,sep=".",remove = FALSE)
all.30.p$trt.load<-as.factor(all.30.p$trt.load)
all.30.p$trt.load<-factor(all.30.p$trt.load,levels=c("0.y","10.y","10.n"))


fd.plot3<-ggplot(all.30.p, aes(x=ttend,group=interaction(temp.var,load.0),fill=trt.load)) 
fd.plot3+geom_density(lwd=1
)+scale_fill_manual(values = alpha(c("#56B4E9","#D55E00","red"), 0.6),
                    name=c("Emergence?"),
                    breaks=c("0.y","10.y","10.n"),
                    labels=c("30+/-0, yes","30+/-10, yes","30+/-10, none")
)+labs(x="Age [days]",y="Density"
)+theme(text = element_text(family=("Cambria")),
        axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        legend.key.width=unit(15,"mm"),
        legend.key.height = unit(10,"mm"),
        legend.text=element_text(size=24),
        legend.title=element_text(size=26),
        legend.position = c(.6,.8),
        legend.background = element_rect(color="black",linetype="solid",size=1))



#---------

#Preparing wasp data for development time to each stage plot

tv.para<-subset(tv.ab,treatment=="para")

tv.para$ovp.day<-1

cc.long<-gather(tv.para,stage,devtime,waspdev.tot,waspdev.int,ovp.day)

cc.long$stage<-gsub("waspdev.tot", "eclosion",cc.long$stage)
cc.long$stage<-gsub("waspdev.int", "emergence",cc.long$stage)
cc.long$stage<-gsub("ovp.day", "oviposition",cc.long$stage)

cc.long$devtime[(cc.long$devtime=="1")]<-0

cc.long$stage<-factor(cc.long$stage,levels=c("oviposition","emergence","eclosion"))


dev.sum<-summarySE(cc.long, measurevar = "devtime", groupvars = c("temp.avg","temp.var","stage"),na.rm=TRUE)
dev.sum



waspdev.plot<-ggplot(dev.sum,aes(x=stage,y=devtime,group=temp.var,color=temp.var))
waspdev.plot+geom_point(size=4,shape=17
)+geom_line(size=1.4
)+geom_errorbar(aes(ymin=devtime-se,ymax=devtime+se),width=.5,size=1
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10")
)+scale_x_discrete(label=c("Egg","Emerge","Eclose")
)+labs(x="Stage",y="Age [days]"
)+facet_wrap(~temp.avg
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=24),
        axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.key.width=unit(15,"mm"),
        legend.background = element_rect(color="black",linetype="solid"))


#---------

#Wasp survivorship curve

tv.para$suc.ovp[is.na(tv.para$suc.ovp)]<-0
tv.para$suc.ovp<-ifelse(tv.para$temp.avg=="30" & tv.para$temp.var=="10",.5,tv.para$suc.ovp)
tv.para<-subset(tv.para,suc.ovp!="0")
tv.para$suc.ovp[(tv.para$suc.ovp==".5")]<-NA


tv.para$stnd.load<-(tv.para$load/tv.para$load)
tv.para$stnd.load[is.nan(tv.para$stnd.load)]<-NA

tv.para$stnd.em<-(tv.para$num.em/tv.para$load)
tv.para$stnd.em[is.nan(tv.para$stnd.em)]<-NA

tv.para$stnd.coc<-(tv.para$num.coc/tv.para$load)
tv.para$stnd.coc[is.nan(tv.para$stnd.coc)]<-NA

tv.para$stnd.ecl<-(tv.para$num.ecl/tv.para$load)
tv.para$stnd.ecl[is.nan(tv.para$stnd.ecl)]<-NA


para.sub<-select(tv.para,bug.id,treatment,temp.avg,temp.var,stnd.load,stnd.em,stnd.coc,stnd.ecl)

para.sub$stnd.load[is.nan(para.sub$stnd.load)]<-0
para.sub$stnd.em[is.nan(para.sub$stnd.em)]<-0
para.sub$stnd.coc[is.nan(para.sub$stnd.coc)]<-0
para.sub$stnd.ecl[is.nan(para.sub$stnd.ecl)]<-0


para.sub<-gather(para.sub,stage,prop.surv,stnd.load,stnd.em,stnd.coc,stnd.ecl)

para.sub$stage<-factor(para.sub$stage,levels=c("stnd.load","stnd.em","stnd.coc","stnd.ecl"))


surv.sum<-summarySE(para.sub,measurevar = "prop.surv",groupvars = c("temp.avg","temp.var","stage"),na.rm=TRUE)
surv.sum

surv.plot2<-ggplot(surv.sum,aes(x=stage,y=prop.surv,group=temp.var,color=temp.var))
surv.plot2+geom_point(size=4,shape=17
)+geom_line(size=1.4
)+geom_errorbar(aes(ymin=prop.surv-se,ymax=prop.surv+se),width=.5,size=1
)+scale_color_manual(values=c("#56B4E9","#D55E00"),name=c("Fluctuation [C]"),
                     breaks=c("0","10"),labels=c("0","10")
)+scale_x_discrete(label=c("Eggs","Emerge","Pupa","Eclose")
)+labs(x="Stage",y="Survival [%]"
)+facet_wrap(~temp.avg
)+theme(text = element_text(family=("Cambria")),
        strip.background = element_rect(colour="black",linetype = "solid",fill="white",
                                        size = 1),
        strip.text = element_text(size=24),
        axis.line.x=element_line(colour = 'black', size = 1),
        axis.line.y=element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = 'black', size = 1),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.key.width=unit(15,"mm"),
        legend.key.height = unit(10,"mm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.position = c(.8,.8),
        legend.background = element_rect(color="black",linetype="solid",size=1))




#-----------------------

#Analysis

#Remove temp.var 5 from the data set. Also removed unsuc ovp

tv.long.no5<-subset(tv.long,temp.var!=5)
tv.long.no5$temp.var<-factor(tv.long.no5$temp.var)


##MODEL 1--no +/- 5 fluctuation
##log.mass is numeric, day.age is integer, treat and temp.avg and temp.var are factors
##random effect is linear day age by individual

#try model selection with factor vs numeric for temperature metrics

lms.mod1<-lme(log.mass~(day.age+I(day.age^2)):(temp.var+treatment+temp.avg)^2+temp.avg,
              random=~day.age|bug.id,
              data=tv.long.no5,
              na.action=na.omit,
              method="ML")

anova(lms.mod1)
summary(lms.mod1)



##subset to only parasitized treatment

tv.para<-subset(tv.ab,treatment=="para")


#Making a column for total died (load-num.ecl)

tv.para$tot.died<-tv.para$load-tv.para$num.ecl


#Making tot.died==NA for mongos (load==0)

tv.para$mongo<-ifelse(tv.para$treatment=="para" & tv.para$date.em.j==0, 1,0)

tv.para$tot.died<-ifelse(tv.para$load=="0", 1.5, tv.para$tot.died)
tv.para$tot.died[tv.para$tot.died=="1.5"]<-NA

tv.para$load[tv.para$load==0]<-NA

tv.para$num.ecl<-ifelse(tv.para$mongo==1, 1.5,tv.para$num.ecl)

tv.para$num.ecl[tv.para$num.ecl=="1.5"]<-NA



#making a glm model of wasp total survival to test for overdispersion
##temp.avg==factor, temp.var==factor, load==numeric
###overdispersion is high, should add random effect of individual

wtots.mod1<-glm(cbind(tot.died,num.ecl)~temp.avg*temp.var*load,
                family=quasibinomial,
                data=tv.para,
                na.action = na.omit)

anova(wtots.mod1,test="F")
summary(wtots.mod1)


#Rescaling load

tv.para$resc.ld<-rescale(tv.para$load,to=c(0,1))

#Making a glmer (binomial) model of wasp total survival
##temp.avg==factor, temp.var==numeric, resc.load==numeric, mongos treated as NAs until I determine a better way to deal with them
###won't run without rescaling load
###as of 2.27.18, won't run even when load is rescaled. Don't know what's going on

wtots.mod2<-glmer(cbind(tot.died,num.ecl)~temp.avg*temp.var+resc.ld+(1|bug.id),
                  family=binomial,
                  data=tv.para,
                  na.action = na.omit,
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

#bootstrap of confindence intervals

anova(wtots.mod2,test="Chisq")
summary(wtots.mod2)



#modelling wasp survival to emergence   This model runs!

tv.para$num.em[tv.para$num.em==0]<-NA
tv.para$num.unem[tv.para$num.unem==0]<-NA

wems.mod1<-glmer(cbind(num.unem,num.em)~temp.avg*temp.var*resc.ld+(1|bug.id),
                 family=binomial,
                 data=tv.para,
                 na.action=na.omit,
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))

anova(wems.mod1)
summary(wems.mod1)





#-------

#Early heat shock stuff, if time

ehs <- read_csv("~Ms+Cc_EHS_incomplete_clean.csv")
View(ehs)

ehs.long <- read_csv("~Ms+Cc_EHS_incomplete_clean_long.csv")
View(ehs.long)



#Looking at mongo.age and mongo.mass at the different heat shock treatments

ehs<-subset(ehs,id!="69")


#Making a column "class"--either em, sm.mongo, mongo or wand

ehs$wand[is.na(ehs$wand)]<-0

ehs$date.em.j[is.na(ehs$date.em.j)]<-0
ehs$date.end.j[is.na(ehs$date.end.j)]<-0

ehs$class<-ifelse(ehs$date.em.j>0,"em",
                  ifelse(ehs$mongo.mass==1, "mongo", "sm.mongo"))

ehs$class2<-ifelse(ehs$wand==1, "wand",0)

ehs$class2[ehs$class2==0]<-NA

ehs$class<-coalesce(ehs$class2, ehs$class)



#Plotting % of treatment that falls into each class--that survived, have removed ind that died

#Dropping individuals with incomplete data (still in expt)

ehs.fin<-subset(ehs, date.3.j<=53)
ehs.fin<-subset(ehs.fin, id!=212)

#Calculating % class in each treatment
tot<-ehs.fin %>% count(hs.temp,hs.num)

n.cl<-ehs.fin %>% count(hs.temp,hs.num,class)

#Adding groups with 0 counts
hs.temp<-c(0,0,40,42,42)
hs.num<-c(0,0,1,1,3)
class<-c("sm.mongo","mongo","mongo","mongo","em")
n<-c(0,0,0,0,0)

add<-data.frame(hs.temp,hs.num,class,n)

n.cl<-rbind(n.cl,add)

#Adding a column with total n for each treatment to n.cl

tot.n<-c(23,23,
         42,42,42,
         44,44,44,44,
         42,42,42,42,
         40,40,40,40,
         42,42,42,
         38,38,38,38,
         45,45,45,
         40,40,40,40,
         23,23,
         42,
         42,
         45)

n.cl$tot.n<-tot.n

#Calculating % of class in each treatment

n.cl<-n.cl %>% mutate(perc.class=n/tot.n)


#Plotting perc.class for each treatment

class.plot<-ggplot(n.cl,aes(x=hs.num,y=perc.class,fill=class))
class.plot+geom_bar(position="fill",stat="identity"
)+facet_wrap(~hs.temp)









