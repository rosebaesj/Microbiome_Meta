#*********************************

library(metafor)
library(meta)
library(dmetar)
library(ggplot2)

dev.new(width=5, height=3, unit="cm")
dev.new(width=3, height=3, unit="cm")
data=read.csv("~/Desktop/Journal of Pain/R/forest plot data subgroups wo chen 복사본.csv")
data2=read.csv("~/Desktop/Journal of Pain/R/MCDR&pain with cro w percent SMD (1).csv")
summary=read.csv("~/Desktop/Journal of Pain/R/summary.csv")
#'#ae5da1','#eb6877','#f19149', '#80c269', '#448aca', '#5f52a0'
#SMD SE

ggplot(data = data) + 
geom_segment(aes(x=c(5.85), y=c(0), xend=c(5.85), yend=c(4)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+#R로 계산했을 때는 5.5941 [4.3587; 6.8296]
geom_segment(aes(x=c(5.85), y=c(0), xend=c(-1), yend=c((5.85+1)/1.96)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(5.85), y=c(0), xend=c(5.85+1.96*4), yend=c(4)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(0), y=c(0), xend=c(0), yend=c(4)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))+
  geom_point(data = data, aes(x = SMD, y = SE, colour=Model),show.legend = FALSE,size=1)+
  scale_color_manual(values=c( '#ae5da1', '#5f52a0')) +  
scale_y_reverse("SE of SMD",breaks=c(0.0, 1.0, 2.0, 3.0, 4.0),  labels = scales::number_format(accuracy = 0.1), expand=c(0,0))	+
scale_x_continuous("SMD", breaks=c(0,5,10,15,20),  limits=c(-1, 20), expand=c(0,0))+
theme_classic()


#SMD inversen
ggplot(data = data) + 
geom_segment(aes(x=c(5.85), y=c(0), xend=c(5.85), yend=c(0.5)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(0), y=c(0), xend=c(0), yend=c(0.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))+
  geom_point(data = data, aes(x = SMD, y = inversen, colour=Model),show.legend = FALSE, size=1)+
  scale_color_manual(values=c( '#ae5da1', '#5f52a0')) +  
scale_y_reverse("1/√ (ESS)", expand=c(0,0))	+
scale_x_continuous("SMD", breaks=c(0,5,10,15,20),  limits=c(-1, 20), expand=c(0,0))+
theme_classic()

#legend
ggplot(data = data) + 
geom_segment(aes(x=c(5.85), y=c(0), xend=c(5.85), yend=c(0.5)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(0), y=c(0), xend=c(0), yend=c(0.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))+
  geom_point(data = data, aes(x = SMD, y = inversen, colour=Model) ,size=1)+
  scale_color_manual(limits=c("Normal", "Model"),values=c( '#5f52a0','#ae5da1')) +  
scale_y_reverse("1/√ (ESS)", expand=c(0,0))	+
scale_x_continuous("SMD", breaks=c(0,5,10,15,20),  limits=c(-1, 20), expand=c(0,0))+
theme_classic()

ggplot(data = data) + 
geom_segment(aes(x=c(5.85), y=c(0), xend=c(5.85), yend=c(0.5)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(0), y=c(0), xend=c(0), yend=c(0.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))+
  geom_point(data = data, aes(x = SMD, y = inversen, colour=Model))+
  scale_color_manual(limits=c("Without Cromolyn", "With Cromolyn"),values=c( '#448aca','#13b5b1')) +  
scale_y_reverse("1/√ (ESS)", expand=c(0,0))	+
scale_x_continuous("SMD", breaks=c(0,5,10,15,20),  limits=c(-1, 20), expand=c(0,0))+
theme_classic()

#regression 돌린것들
ggplot(data = data2) + 
geom_segment(aes(x=c(2.28), y=c(0), xend=c(2.28), yend=c(1.5)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(2.28), y=c(0), xend=c(2.28-1.96*1.5), yend=c(1.5)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(2.28), y=c(0), xend=c(2.28 +1.96*1.5), yend=c(1.5)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(0), y=c(0), xend=c(0), yend=c(1.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))+
  geom_point(data = data2, aes(x = SMD, y = SE, colour=Cromolyn),show.legend = FALSE,size=1)+
  scale_color_manual(values=c('#13b5b1','#448aca')) +  
scale_y_reverse("SE of SMD", expand=c(0,0))	+
scale_x_continuous("SMD", limits=c(-2, 6), expand=c(0,0))+
theme_classic()


ggplot(data = data2) + 
geom_segment(aes(x=c(2.28), y=c(0), xend=c(2.28), yend=c(0.5)), linetype=2, size=0.5, colour="grey70", expand=c(0,0))+
geom_segment(aes(x=c(0), y=c(0), xend=c(0), yend=c(0.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))+
  geom_point(data = data2, aes(x = SMD, y = inversen, colour=Cromolyn),show.legend = FALSE,size=1)+
  scale_color_manual(values=c('#13b5b1','#448aca')) +  
scale_y_reverse("1/√ (ESS)", expand=c(0,0))	+
scale_x_continuous("SMD",  limits=c(-2, 6), breaks=c(-2, 0,2, 4, 6), expand=c(0,0))+
theme_classic()


 dev.new(width=3.5, height=3, unit="cm")
#Study 별 *********************************
ggplot(data=data2, aes(x=MCDR, y=SMD, weight=weight))+
	geom_smooth(method = lm, stat = "smooth", color='black', fill="grey70",size=0.5, level=0.95, formula=y~x) +
   geom_point(shape =1, aes(size=ESS, color=Study, stroke=0.7), show.legend = FALSE) + #show.legend = FALSE
   scale_color_manual(values=c('#ae5da1','#eb6877','#f19149', '#80c269', '#448aca', '#5f52a0')) + # 원하는 색상 코드로 변경 
   theme_classic() +
    scale_size(range=c(5,10), breaks=c(5,10), labels=c(5, 10))+
scale_x_continuous("Mast Cell Degranulation Ratio (%)", breaks=c(20, 40, 60, 80),  limits=c(20, 80), expand=c(0,0))+
 scale_y_continuous("Pain Threshoold (SMD)", breaks=c(-1,0,1,2,3,4,5), limits=c(-1, 5.5), expand=c(0,0))



 #Group별*********************************
 ggplot(data=data2, aes(x=MCDR, y=SMD, weight=weight))+
	geom_smooth(method = lm, stat = "smooth", color='black', fill="grey70",size=0.5, level=0.95, formula=y~x) +
   geom_point(aes(size=ESS, colour=Subgroup,alpha=0.99), show.legend = FALSE) +
   scale_color_manual(limits=c("IA", "MA", "MA+cromolyn", "EA", "EA+cromolyn"), values=c('#ae5da1', '#eb6877','#f19149', '#448aca', '#80c269' )) + # 원하는 색상 코드로 변경 
       theme_classic() +
       scale_size(range=c(5,10), breaks=c(5,10), labels=c(5, 10)) +
scale_x_continuous("Mast Cell Degranulation Ratio (%)", breaks=c(20, 40, 60, 80),  limits=c(20, 80), expand=c(0,0))+
 scale_y_continuous("Pain Threshoold (SMD)", breaks=c(-1,0,1,2,3,4,5), limits=c(-1, 5.5), expand=c(0,0))

#!!!!!regression tests
data=read.csv("~/Desktop/Journal of Pain/R/forest plot data subgroups wo chen 복사본.csv")
data2=read.csv("~/Desktop/Journal of Pain/R/MCDR&pain with cro w percent SMD (1).csv")

 y.se=rma(SMD,SE,data=data)
 y.n=rma(SMD,inversen,data=data)
 y2.se=rma(SMD,SE,data=data2)
 y2.n=rma(SMD,inversen,data=data2)

 regtest(y.se,y.se$SMD,y.se$SE)
 regtest(y.n,y.n$SMD,y.n$inversen)
 regtest(y2.se,y2.se$SMD,y2.se$SE)
 regtest(y2.n,y2.n$SMD,y2.n$inversen)

m.pain=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Study,data=data2,byvar=MCDR,sm="SMD")
r.pain=metareg(m.pain)

##funnel plot drawing
d=read.csv("~/Desktop/Journal of Pain/R/forest plot data.csv")
y=rma(MD,SEMD,data=d)
funnel(y,steps=4,ylim=c(3,0))

d=read.csv("~/Desktop/Journal of Pain/R/forest plot data subgroups wo chen.csv")
y=rma(MD,inversen,data=d)
funnel(y,steps=4,ylim=c(1,0))
	

d=read.csv("~/Desktop/Journal of Pain/R/forest plot data subgroups.csv")
y=rma(MD,size,data=d)
funnel(y)
funnel(y, ylim=c(1,0))





#********************
estimate = 5.85
se = 0.67857143

se.seq=seq(0, max(dat$corr_zi_se), 0.001)

ggplot(data = d) + 

  geom_point(data = d, aes(x = Size, y = Mean), size = 2,
        colour = "black", shape = 21,fill = filling3) + 
  ylim(0, 8)
  geom_line(data = LINE3, aes(x = 1:(max(CG_PLOT1$Size) + 25), 
        y = M3 + qnorm(0.975) * SD3 / N3), size = 1, colour = "grey70",
        linetype = 5) +
  geom_line(data = LINE3, aes(x = 1:(max(CG_PLOT1$Size) + 25), 
        y = M3 - qnorm(0.975) * SD3 / N3), size = 1, colour = "grey70",
        linetype = 5) +
  geom_segment(xend = max(CG_PLOT1$Size)+25,yend=mean(LINE3$M3,na.rm=T)),
       aes(x = 1, y = mean(LINE3$M3,na.rm=T), size=1, colour="grey70") +
theme_classic() +


dev.new(width=5, height=3, unit="cm")


#MD SEMD
ggplot(data = data) + 
  geom_point(data = data, aes(x = MD, y = SEMD, colour=Model,))+
  scale_color_manual(values=c('#80c269', '#eb6877','#f19149',  '#448aca', '#ae5da1')) +  
scale_y_reverse("SEMD", expand=c(0,0))	+
scale_x_continuous(limits=c(0, 80), expand=c(0,0))+
theme_classic()+
geom_segment(aes(x=c(34.44), y=c(0), xend=c(34.44), yend=c(5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))+# R로 계산시 34.3999 [25.4545; 43.3452] 

#MD inversen
ggplot(data = data)
  geom_point(data = data, aes(x = MD, y = inversen, colour=Model,))+
  scale_color_manual(values=c('#80c269', '#eb6877','#f19149',  '#448aca', '#ae5da1')) +  
scale_y_reverse("inverse sqare root of n", expand=c(0,0))	+
scale_x_continuous(limits=c(0, 80), expand=c(0,0))+
theme_classic()+
geom_segment(aes(x=c(34.44), y=c(0), xend=c(34.44), yend=c(0.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))

 
  
  
  
  
  
  
  
  
  data=read.csv("~/Desktop/Journal of Pain/R/.csv")
    
#regression 돌린것****************
ggplot(data = data) + 
  geom_point(data = data, aes(x = SMD, y = SE, colour=Cromolyn,))+
  scale_color_manual(values=c('#80c269', '#eb6877','#f19149',  '#448aca', '#ae5da1')) +  
scale_y_reverse("SE of SMD", expand=c(0,0))	+
scale_x_continuous(breaks=c(-5,0,5,10),  limits=c(-5, 10), expand=c(0,0))+
theme_classic()+
geom_segment(aes(x=c(2.28), y=c(0), xend=c(2.28), yend=c(1.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))#R로 계산했을 때는 5.5941 [4.3587; 6.8296]

ggplot(data = data) + 
  geom_point(data = data, aes(x = SMD, y = inversen, colour=Cromolyn,))+
  scale_color_manual(values=c('#80c269', '#eb6877','#f19149',  '#448aca', '#ae5da1')) +  
scale_y_reverse("inversen", expand=c(0,0))	+
scale_x_continuous(breaks=c(-5,0,5,10),  limits=c(-5, 10), expand=c(0,0))+
theme_classic()+
geom_segment(aes(x=c(2.28), y=c(0), xend=c(2.28), yend=c(0.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))

ggplot(data = data) + 
  geom_point(data = data, aes(x = MD, y = SEMD, colour=Cromolyn,))+
  scale_color_manual(values=c('#80c269', '#eb6877','#f19149',  '#448aca', '#ae5da1')) +  
scale_y_reverse("SEMD", expand=c(0,0))	+
scale_x_continuous(breaks=c(-5,0,5,10),  limits=c(-5, 10), expand=c(0,0))+
theme_classic()+
geom_segment(aes(x=c(2.26), y=c(0), xend=c(2.26), yend=c(5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))# R로 계산시 34.3999 [25.4545; 43.3452] 

ggplot(data = data) + 
  geom_point(data = data, aes(x = MD, y = inversen, colour=Cromolyn,))+
  scale_color_manual(values=c('#80c269', '#eb6877','#f19149',  '#448aca', '#ae5da1')) +  
scale_y_reverse("inverse sqare root of n", expand=c(0,0))	+
scale_x_continuous(breaks=c(-5,0,5,10),  limits=c(-5, 10), expand=c(0,0))+
theme_classic()+
geom_segment(aes(x=c(2.26), y=c(0), xend=c(2.26), yend=c(0.5)), linetype=1, size=0.5, colour="grey70", expand=c(0,0))

panel.grid) 

#Study 별 *********************************
ggplot(data, aes(x=MCDR, y=SMD, weight=size))+
	geom_smooth(method = lm, stat = "smooth", color='black', fill="#e5e5e5",size=0.5, level=0.95, formula=y~x) +
   geom_point(shape =1, aes(size=size, color=Study, stroke=0.7)) +
   scale_color_manual(values=c('#ae5da1','#eb6877','#f19149', '#80c269', '#448aca', '#5f52a0')) + # 원하는 색상 코드로 변경 
   theme_classic() +
    scale_size(range=c(5,10), breaks=c(5,10), labels=c(5, 10))+
scale_x_continuous("Mast Cell Degranulation Ratio (%)", breaks=c(20, 40, 60, 80),  limits=c(20, 80), expand=c(0,0))+
 scale_y_continuous("Pain Threshoold (SMD)", breaks=c(-1,0,1,2,3,4,5), limits=c(-1, 5.5), expand=c(0,0))+
 dev.new(width=5.5, height=3, unit="cm")

 
 
 
 #Group별*********************************
 ggplot(data, aes(x=MCDR, y=SMD))+
	geom_smooth(method = lm, stat = "smooth", color='black', fill="#e5e5e5",size=0.5, level=0.95, formula=y~x) +
   geom_point(aes(size=size, colour=Subgroup,alpha=0.95)) +
   scale_color_manual(values=c('#eb6877','#f19149',  '#80c269', '#448aca', '#ae5da1')) + # 원하는 색상 코드로 변경 
       theme_classic() +
       scale_size(range=c(5,10), breaks=c(5,10), labels=c(5, 10)) +
scale_x_continuous("Mast Cell Degranulation Ratio (%)", breaks=c(20, 40, 60, 80),  limits=c(20, 80), expand=c(0,0))+
 scale_y_continuous("Pain Threshoold (SMD)", breaks=c(-1,0,1,2,3,4,5), limits=c(-1, 5.5), expand=c(0,0))
 dev.new(width=5.5, height=3, unit="cm")


#SMD, SE of SMD***************
d=read.csv("~/Desktop/R/forest plot data subgroups wo chen.csv")
y=rma(SMD,SE,data=d)
dev.new(width=5, height=4, unit="cm")
funnel(y,steps=4,ylim=c(3,0), xlab="SMD", ylab="SE of SMD")

#MD, SE of MD
d=read.csv("~/Desktop/R/forest plot data subgroups wo chen.csv")
y=rma(MD,SEMD,data=d)
dev.new(width=5, height=4, unit="cm")
funnel(y,steps=4,ylim=c(3,0), xlab="MD", ylab="SE of MD")
funnel(y,steps=5,ylim=c(1,0), xlab="MD", ylab="inverse square root of n", back="white", hlines="grey90", refline=c(2.98,8.76))


#SMD, inverse square root of n
d=read.csv("~/Desktop/R/forest plot data subgroups wo chen.csv")
y=rma(SMD,inversen,data=d)
dev.new(width=5, height=4, unit="cm")
funnel(y,steps=5,ylim=c(1,0), xlab="SMD", ylab="inverse square root of n")

#MD, inverse square root of n
d=read.csv("~/Desktop/R/forest plot data subgroups wo chen.csv")
y=rma(MD,inversen,data=d)
dev.new(width=5, height=4, unit="cm")
funnel(y,steps=5,ylim=c(1,0), xlab="MD", ylab="inverse square root of n", level=0, back="white", hlines="grey90", col=(d$col2))



nd=read.csv("~/Desktop/R/forest plot data normal.csv")
ny=rma(MD,SEMD,data=nd)
funnel(ny,steps=4,ylim=c(3,0))

md=read.csv("~/Desktop/R/forest plot data model.csv")
my=rma(MD,SEMD,data=md)
funnel(my,steps=4,ylim=c(3,0))

##regression test, rank test
regtest(y,y$MD,y$SEMD)
ranktest(y)

#metacont 기본 데이터, Fixed 인지 Random인지, heterogeneity Q 값도 나옴 보려면 a 만 치면 됨

a=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Study.or.Subgroup,data=d)
na=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Study.or.Subgroup,data=nd)
ma=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Study.or.Subgroup,data=md)

#forestplot
forest(a)
forest(na)
forest(ma)

#forestplot 디자인
forest(a, col.study = "black", col.square = "blue",  col.inside = "white", col.diamond = "gray",  col.diamond.fixed = col.diamond,
  col.diamond.random = col.diamond,
  col.diamond.lines = "black",
  col.diamond.lines.fixed = col.diamond.lines,
  col.diamond.lines.random = col.diamond.lines,  print.Q = FALSE)

#bubbleplot
library(metafor)
library(meta)
pain=read.csv("~/Desktop/R/MCDR&pain with cro.csv")
m.pain=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Study,data=pain,byvar=MCDR,sm="SMD")
r.pain=metareg(m.pain)
bubble(r.pain, xlab = "Publication Year",col.line = "blue",studlab = TRUE)
bubble(r.pain,col.line = "blue")

#subgroup analysis
library(dmetar)
d=read.csv("~/Desktop/R/forest plot data subgroups.csv")
a=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Study.or.Subgroup,data=d,sm="SMD")
subgroup.analysis.mixed.effects(x = a, subgroups = d$Model)
subgroup.analysis.mixed.effects(x = a, subgroups = d$by2010)
subgroup.analysis.mixed.effects(x = a, subgroups = d$Acupuncture)
subgroup.analysis.mixed.effects(x = a, subgroups = d$Animal)

#SMD
a=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Study.or.Subgroup,data=d, method.smd = gs("method.smd"))
subgroup.analysis.mixed.effects(x = a, subgroups = d$Model)
subgroup.analysis.mixed.effects(x = a, subgroups = d$by2010)
subgroup.analysis.mixed.effects(x = a, subgroups = d$Acupuncture)
subgroup.analysis.mixed.effects(x = a, subgroups = d$Animal)

#year
a=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Study.or.Subgroup,data=d,byvar=Y)
r=metareg(a)
bubble(r, xlab = "Publication Year",col.line = "blue",studlab = TRUE)
bubble(rn,col.line = "blue")


# 데이터 불러오기 
data <- read.csv("~/Desktop/R/MCDR&pain with cro w percent SMD.csv", header=TRUE)

# 라이브러리 호출
install.packages("ggplot2") #최초 1회만 하면 됨 
library(ggplot2)

 
 #Study 별 *********************************
ggplot(data, aes(x=MCDR, y=SMD))+
	geom_smooth(method = lm, stat = "smooth", color='black', fill="#e5e5e5",size=0.5, level=0.95, formula=y~x) +
   geom_point(shape =1, aes(size=size, color=Study, stroke=0.7)) +
   scale_color_manual(values=c('#ae5da1','#eb6877','#f19149', '#80c269', '#448aca', '#5f52a0')) + # 원하는 색상 코드로 변경 
   theme_classic() +
    scale_size(range=c(5,10), breaks=c(5,10), labels=c(5, 10))+
scale_x_continuous("Mast Cell Degranulation Ratio (%)", breaks=c(20, 40, 60, 80),  limits=c(20, 80), expand=c(0,0))+
 scale_y_continuous("Pain Threshoold (SMD)", breaks=c(-1,0,1,2,3,4,5), limits=c(-1, 5.5), expand=c(0,0))+
 dev.new(width=5.5, height=3, unit="cm")

 
 
 
 #Group별*********************************
 ggplot(data, aes(x=MCDR, y=SMD))+
	geom_smooth(method = lm, stat = "smooth", color='black', fill="#e5e5e5",size=0.5, level=0.95, formula=y~x) +
   geom_point(aes(size=size, colour=Subgroup,alpha=0.95)) +
   scale_color_manual(values=c('#eb6877','#f19149',  '#80c269', '#448aca', '#ae5da1')) + # 원하는 색상 코드로 변경 
       theme_classic() +
       scale_size(range=c(5,10), breaks=c(5,10), labels=c(5, 10)) +
scale_x_continuous("Mast Cell Degranulation Ratio (%)", breaks=c(20, 40, 60, 80),  limits=c(20, 80), expand=c(0,0))+
 scale_y_continuous("Pain Threshoold (SMD)", breaks=c(-1,0,1,2,3,4,5), limits=c(-1, 5.5), expand=c(0,0))
 dev.new(width=5.5, height=3, unit="cm")
 
 #*********************************
 a=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Subgroup,data=data,byvar=Y,sm="SMD")
r=metareg(a)
 
 
  ggplot(data, aes(x=MCDR, y=SMD))+
ylim(-1,5.5)+
	geom_smooth(method = lm, stat = "smooth", color='black', fill="#e5e5e5",size=0.5, level=0.95, formula=y~x) +
   geom_point(aes(size=size, colour=Study, alpha=0.95)) +
   scale_color_manual(values=c('#eb6877','#f19149',  '#80c269', '#448aca', '#ae5da1')) + # 원하는 색상 코드로 변경 
   scale_size(range=c(10,20)) +
   theme_classic() +
   xlab("Mast Cell Degranulation Ratio (%)") + #x축 이름 
   ylab("Pain Threshoold(Standard Mean Difference)") #y축 이름 
 
 
 
 
 
 
 
 
 
 #MD 
 data <- read.csv("~/Desktop/R/MCDR&pain with cro MD.csv", header=TRUE)
 ggplot(data, aes(x=MCDR, y=MD))+
	geom_smooth(method = lm, stat = "smooth", color='black', fill="#e5e5e5",size=0.5, level=0.95, formula=y~x) + #회귀식 라인
   geom_point(aes(size=size, colour=Study), alpha=0.7) +
   scale_color_manual(values=c('#f29c9f','#f6b37f', '#fff799', '#acd598', '#7ecef4', '#8f82bc')) + # 원하는 색상 코드로 변경 
   scale_size(range=c(5,20)) +
   theme_classic() +
   xlab("Mast Cell Degranulation Ratio (%)") + #x축 이름 
   ylab("Pain Threshoold(Mean Difference)") + #y축 이름 
   geom_smooth(method = lm, stat = "smooth", color='black', fill="#e5e5e5",size=0.5, level=0.95, formula=y~x) #회귀식 라인

  #year
a=metacont(aTotal,amean,aSD,cTotal,cmean,cSD,Subgroup,data=data,byvar=Y,sm="SMD")
r=metareg(a)
bubble(r, xlab = "Publication Year",col.line = "blue",studlab = TRUE)
bubble(rn,col.line = "blue")
 


















