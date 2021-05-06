rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(colorscience)
library(scatterplot3d)
library(rgl)
library(geometry)
library(matlib)
library(tikzDevice)
library(readr)

data_dir<-"./"

#Stiles & Burch (1955) 2-deg, RGB CMFs
sb.df<-read.csv(paste0(data_dir,"sbrgb2.csv"),header=F)
colnames(sb.df)<-c("wavelength","r.lambda","g.lambda","b.lambda")
sb.df
#CIE 1931 2-deg, XYZ CMFs
xyz.df<-read.csv(paste0(data_dir,"ciexyz31.csv"),header=F)
colnames(xyz.df)<-c("wavelength","x.lambda","y.lambda","z.lambda")
xyz.df
#2-deg XYZ CMFs transformed from the CIE (2006) 2-deg LMS cone fundamentals
o.df<-read.csv(paste0(data_dir,"linss2_10e_5.csv"),header=F)
colnames(o.df)<-c("wavelength","l.lambda","m.lambda","s.lambda")
o.df

o.df<-o.df %>% inner_join(sb.df,by="wavelength") %>% inner_join(xyz.df,by="wavelength")
o.df

pho.df<-read.csv(paste0(data_dir,"vljve.csv"),header=F)
sco.df<-read.csv(paste0(data_dir,"scvle.csv"),header=F)

colnames(pho.df)<-c("wavelength","V")
colnames(sco.df)<-c("wavelength","Vp")

purple1.wl<-400
purple2.wl<-680
lms.df <- o.df %>% filter(wavelength>=purple1.wl,wavelength<=purple2.wl )
lms.df

lms.purple1<-lms.df %>% filter(wavelength==purple1.wl)
lms.purple2<-lms.df %>% filter(wavelength==purple2.wl)

purple.funcs<-apply(rbind(lms.purple1,lms.purple2), 2, approxfun, x=c(purple1.wl, purple2.wl))
purples.lms.df<-data.frame(t(sapply(seq(purple1.wl,purple2.wl,30),function(x){c(NA,purple.funcs$l.lambda(x),purple.funcs$m.lambda(x),purple.funcs$s.lambda(x))})))
colnames(purples.lms.df)=colnames(lms.df)
lms.df<-rbind(lms.df,purples.lms.df %>% arrange(-row_number()))

lms.df

lmss<-as.list(data.frame(apply(lms.df[,c("l.lambda","m.lambda","s.lambda")],1,c)))
names(lmss)<-lms.df$wavelength
lmss


#http://www.cvrl.org/database/text/cienewxyz/cie2012xyz2.htm
lmst<-matrix(c( 1.94735469,  -1.41445123,   0.36476327,
                0.68990272,   0.34832189,            0,
                0,            0,            1.93485343), byrow=T,3,3)


#https://en.wikipedia.org/wiki/SRGB
xyzrgbt<-t(matrix(c( 3.2406, -1.5372, -0.4986, 
                     -0.9689,  1.8758,  0.0415,
                     0.0557, -0.2040,  1.0570),3,3))

gammacor<-function(x){
  if(x<=0.0031308){
    12.92*x
  } else {
    a<-0.055
    (1+a)*x^(1/2.4)-a
  }
}


lms2xyz<-function(x){lmst %*% x}

lms2rgbtt<-function(x){((xyzrgbt %*% lmst) %*% x)}
xyz.df<-data.frame(do.call("rbind",lapply(lmss,function(x){t(lms2xyz(x))})))
colnames(xyz.df)<-c('X','Y','Z')

lab.df<-t(apply(xyz.df,1,XYZ2Lab))
colnames(lab.df)<-c('L*','a*','b*')
lab.df

rgb.df<-do.call("rbind",lapply(lmss,function(x){pmin(1,pmax(0,sapply(lms2rgbtt(x),gammacor)))}))
colnames(rgb.df)<-c('R','G','B')

spectrum.df<-cbind(lms.df,xyz.df,rgb.df)

lmsnorm.df<-spectrum.df %>% select(Ln=l.lambda,Mn=m.lambda,Sn=s.lambda)
chromaticity.df<-spectrum.df %>% select(x=X,y=Y,z=Z)

lmsnorm.df<-t(apply(lmsnorm.df,1,function(x){x/sum(x)}))
chromaticity.df<-t(apply(chromaticity.df,1,function(x){x/sum(x)}))

spectrum.df<-cbind(spectrum.df,lmsnorm.df,chromaticity.df,lab.df) %>% mutate(wl.labels=as.character(wavelength))
spectrum.df$wl.labels[as.logical(spectrum.df$wavelength %% 15)]=""
spectrum.df$wl.labels[as.logical(spectrum.df$wavelength < 400 | spectrum.df$wavelength > 630 | is.na(spectrum.df$wavelength))]=""

spectrum.df

#this section needs to be run to translate rgl graphs with mouse, need to run this after you create the rgl window
pan3d <- function(button) {
start <- list()

begin <- function(x, y) {
 start$userMatrix <<- par3d("userMatrix")
  start$viewport <<- par3d("viewport")
  start$scale <<- par3d("scale")
   start$projection <<- rgl.projection()
  start$pos <<- rgl.window2user( x/start$viewport[3], 1 - y/start$viewport[4], 0.5, 
                              projection=start$projection)
}

update <- function(x, y) {
 xlat <- (rgl.window2user( x/start$viewport[3], 1 - y/start$viewport[4], 0.5,
                         projection = start$projection) - start$pos)*start$scale
mouseMatrix <- translationMatrix(xlat[1], xlat[2], xlat[3])
par3d(userMatrix = start$userMatrix %*% t(mouseMatrix) )
}
rgl.setMouseCallbacks(button, begin, update)
cat("Callbacks set on button", button, "of rgl device",rgl.cur(),"\n")
}
pan3d(2)

#for grabbing rgl view matrix
zoom<-par3d()$zoom
userMatrix<-par3d()$userMatrix
windowRect<-par3d()$windowRect
zoom
paste(userMatrix,collapse=",")
windowRect


#luminous efficiency
leff.df<-inner_join(sco.df,pho.df,by="wavelength")
leff.df<-inner_join(leff.df,spectrum.df %>% select(c(wavelength,L,M,S,R,G,B)))  %>% gather(cone.type, response, c(V,Vp))
leff.df
ggplot(leff.df, aes(x=wavelength,y=response,group=cone.type))+
  geom_point(size=2,aes(color=rgb(R,G,B)),data=. %>% filter(cone.type=="V"))+
  geom_line(aes(color=rgb(R,G,B)), data=. %>% filter(cone.type=="V"))+
  scale_color_identity()+
  
  geom_point(size=2,aes(color="grey80"),data=. %>% filter(cone.type=="Vp"))+
  geom_line(aes(color="grey80"), data=. %>% filter(cone.type=="Vp"))+
  
  geom_text(data=data.frame(cone.type=c('scoptopic\nluminous\nefficiency','photopic\nluminous\nefficiency'),wavelength=c(400,605),response=c(0.5,0.5)), aes(label=cone.type),size=5,hjust=0,vjust=-1,color="white")+
  labs(x="wavelength (nm)", y="normalized cone sensitivity")+
  theme(
    axis.text = element_text(color="black",size = 12),
    axis.title=element_text(face="bold", size=14),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.minor = element_line(colour = "grey40"),
    panel.background = element_rect(fill = "black")
  )+
  scale_x_continuous(minor_breaks=seq(400,680,25), breaks=seq(400,680,50))



#CMFs

t.df <- spectrum.df %>% select(c(wavelength,r.lambda,g.lambda,b.lambda,x.lambda,y.lambda,z.lambda,l.lambda,m.lambda,s.lambda,R,G,B)) %>% gather(cmf, response, c(c(r.lambda,g.lambda,b.lambda,,x.lambda,y.lambda,z.lambda,l.lambda,m.lambda,s.lambda)))
t.df <- %>% mutate(cmf = case_when(
  startsWith(cmf,"r") ~ "rgb",
  TRUE ~ "other"))

t.df <- t.df %>% mutate(cmf.family = case_when(
  startsWith(cmf,"r") ~ "CIE RGB",
  startsWith(cmf,"g") ~ "CIE RGB",
  startsWith(cmf,"b") ~ "CIE RGB",
  startsWith(cmf,"x") ~ "CIE XYZ",
  startsWith(cmf,"y") ~ "CIE XYZ",
  startsWith(cmf,"z") ~ "CIE XYZ",
  startsWith(cmf,"l") ~ "LMS",
  startsWith(cmf,"m") ~ "LMS",
  startsWith(cmf,"s") ~ "LMS",
  TRUE ~ "other"))

cmf_labels.df<-read_delim(paste0(data_dir,"cmflabels.csv"), delim=',', escape_double=FALSE, escape_backslash=TRUE)

tikzDevice::tikz(file = paste0(data_dir,"cmfs.tex"), width = 5, height = 8,standAlone = T)
ggplot(t.df, aes(x=wavelength,y=response,group=cmf,color=rgb(R,G,B)))+
  facet_grid(rows = vars(cmf.family)) +
  geom_point(size=2)+
  geom_line()+
  scale_color_identity()+
  #geom_text(data=data.frame(cmf=c"L,M,S",wavelength=c(605,505,460),response=c(0.8,0.8,0.8)), aes(label=cmf),size=8,hjust=0,vjust=-1,color="white")+
  labs(x="wavelength (nm)", y="relative response")+
  theme(
    axis.text = element_text(color="black",size = 12),
    axis.title=element_text(face="bold", size=14),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.minor = element_line(colour = "grey40"),
    panel.background = element_rect(fill = "black")
  )+
  scale_x_continuous(minor_breaks=seq(400,680,25), breaks=seq(400,680,50))+
  geom_text(data=cmf_labels.df,aes(label = label,group="none"),size=4,hjust=0,vjust=-1,color="white")
dev.off()


#cone spectral sensitivities
t.df <- spectrum.df %>% select(c(wavelength,l.lambda,m.lambda,s.lambda,R,G,B)) %>% gather(cone.type, response, c(l.lambda,m.lambda,s.lambda))
t.df

ggplot(t.df, aes(x=wavelength,y=response,group=cone.type,color=rgb(R,G,B)))+
geom_point(size=2)+
geom_line()+
scale_color_identity()+
geom_text(data=data.frame(cone.type=c('L','M','S'),wavelength=c(605,505,460),response=c(0.8,0.8,0.8)), aes(label=cone.type),size=8,hjust=0,vjust=-1,color="white")+
labs(x="wavelength (nm)", y="normalized cone sensitivity")+
theme(
    axis.text = element_text(color="black",size = 12),
    axis.title=element_text(face="bold", size=14),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.minor = element_line(colour = "grey40"),
    panel.background = element_rect(fill = "black")
  )+
scale_x_continuous(minor_breaks=seq(400,680,25), breaks=seq(400,680,50))



ggsave("cones.png")



#Lab spectral sensitivities
t.df <- spectrum.df %>% select(c("wavelength",'L*','a*','b*')) %>% gather(primary, response, c('L*','a*','b*'))
t.df

ggplot(t.df, aes(x=wavelength,y=response,group=primary,color=primary))+
  geom_point(size=2)+
  geom_line()
    
    
    geom_text(data=data.frame(cone.type=c('L*','a*','b*'), wavelength=c(605,505,460),response=c(0.8,0.8,0.8)), aes(label=cone.type),size=8,hjust=0,vjust=-1,color="white")+
  labs(x="wavelength (nm)", y="normalized cone sensitivity")+
  theme(
    axis.text = element_text(color="black",size = 12),
    axis.title=element_text(face="bold", size=14),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.minor = element_line(colour = "grey40"),
    panel.background = element_rect(fill = "black")
  )+
  scale_x_continuous(minor_breaks=seq(400,680,25), breaks=seq(400,680,50))


ggsave("lab.png")



#run this 'with' block for the 3D LMS chart
with(spectrum.df, {
 	if( rgl.cur() == 0 ) {
    rgl.open()
    par3d(zoom= 0.6536796, userMatrix=matrix(c(0.972121298313141,0.171159878373146,0.160263866186142,0,0,0.683490395545959,-0.729959487915039,0,-0.23447859287262,0.709609150886536,0.664435565471649,0,0.207534619480768,0.234377896573699,6.63464732111873e-08,1),4,4),windowRect = c( 1539, 102, 3414, 1413 ) )
    rgl.bg(color = "black" )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
 # rgl.viewpoint(theta = 0, phi = -130, zoom = 1)
  axes.title.size=3
  axes3d( edges=c("x--", "y+-", "z--"), col="white", cex=2,nticks=5)
	
	box.color=rgb(0.6,0.6,0.6)
	box.alpha=0.5
	outset=rep(-0.01,4)

    rgl.lines(c(0,0),c(0,1),c(0,0),color="white")
    rgl.lines(c(0,0),c(1,1),c(0,1),color="white")
    rgl.lines(c(0,1),c(1,1),c(0,0),color="white")
    	
    
    rgl.texts(rbind(c(0.5,0,0),c(0,0,0.5)), text = c("L", "S"), color = "white",
             adj = c(6, 6), cex = axes.title.size)
    rgl.texts(rbind(c(1,0.5,0)), text = c("M"), color = "white",
             adj = c(-5, 0), cex = axes.title.size)        
    rgl.points(L,M,S,col=rgb(spectrum.df %>% select(R,G,B)),size=17) 
    text3d(L,M,S,wl.labels,adj=c(-0.3,0),depth_test="always", cex=2)
    
    g.lms.df<-spectrum.df %>% select(L,M,S)
    lmsvolume<-convhulln(g.lms.df)
    triangles3d(g.lms.df[t(lmsvolume), ], color = rgb(1,1,1), depth_test="always",alpha = 0.4)
})

rgl.snapshot("LMS.png")


#For 2D LMS mixing triangle
euler.angles.matrix<-function(phi,theta,psi) {
  D<-t(matrix(c(          cos(phi),      sin(phi),          0,
                         -sin(phi),      cos(phi),          0,
                                 0,              0,         1),3,3))

  C<-t(matrix(c(                 1,             0,           0,
                                 0,    cos(theta),  sin(theta),
                                 0,   -sin(theta),  cos(theta)),3,3))

  B<-t(matrix(c(          cos(psi),  sin(psi),   0,
                         -sin(psi),  cos(psi),   0,
                                 0,         0,   1),3,3))
                                 
  D %*% C %*% B                
}


eangles=c(pi/4,0,0)
basis2d.df<-data.frame(t(apply(matrix(c(1,0,0,
                                        0,1,0,
                                        0,0,1),3,3),
                                        1,function(x){do.call("euler.angles.matrix",as.list(eangles))%*%x})))
colnames(basis2d.df)<-c('L2d','M2d','S2d')
basis2d.df$type=c('L','M','S')
basis2d.df

lms2d.df<-data.frame(t(apply(lmsnorm.df,1,function(x){do.call("euler.angles.matrix",as.list(eangles))%*%x})))
colnames(lms2d.df)<-c('L2d','M2d','S2d')
lms2d.df


lms2rgbttYadj<-function(x){ret<-(lmst %*% x) *c(1,0.5,1)
	                         xyzrgbt %*% ret}


rgb2d.df<-do.call("rbind",lapply(lmss,function(x){pmin(1,pmax(0,sapply(lms2rgbttYadj(x),gammacor)))}))
colnames(rgb2d.df)<-c('R','G','B')
rgb2d.df

spectrum2d.df<-cbind(lms2d.df,rgb2d.df)

ggplot(spectrum2d.df,aes(x=M2d,y=S2d))+
geom_polygon(data=basis2d.df,color="grey50",fill="black",alpha=1)+
geom_point(data=basis2d.df, aes(color="white",size=4))+
geom_text(data=basis2d.df, aes(label=type)) +
geom_polygon(color="grey60",fill="grey60")+
geom_point(aes(color=rgb(R,G,B),size=4))+
scale_color_identity()+
theme(
#    axis.text = element_text(color="black",size = 12),
#    axis.title=element_text(face="bold", size=14),
    axis.title=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "black"),
    legend.position="none"
  )
ggsave("LMStriangle.png")


#run this 'with' block for the 3D XYZ chart
with(spectrum.df, {
 	if( rgl.cur() == 0 ) {
    rgl.open()
    par3d(zoom= 0.5597312, userMatrix=matrix(c(0.85748153924942,0.32283341884613,0.400629490613937,0,0,0.778655052185059,-0.627452194690704,0,-0.514514744281769,0.538028657436371,0.667682349681854,0,-0.00482579399383098,0.195442820474495,3.23043413130719e-08,1),4,4),windowRect = c( 1525, 44, 3525, 1440 ) )
    rgl.bg(color = "black" )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
 	
  box.color=rgb(0.6,0.6,0.6)
  box.alpha=0.5
  outset=rep(-0.01,4)

  rgl.lines(c(0,0),c(0,1),c(0,0),color="white")
  rgl.lines(c(0,0),c(1,1),c(0,1),color="white")
  rgl.lines(c(0,1),c(1,1),c(0,0),color="white")
    	
  axes.title.size=3
  rgl.texts(rbind(c(0.5,0,0),c(0,0,0.5)), text = c("X", "Z"), color = "white",
             adj = c(5, 5), cex = axes.title.size)
  rgl.texts(rbind(c(1,0.5,0)), text = c("Y"), color = "white",
             adj = c(-9, 0), cex = axes.title.size)        
  rgl.points(X,Y,Z,col=rgb(spectrum.df %>% select(R,G,B)),size=17,depth_test="always") 
  text3d(X,Y,Z,wl.labels,adj=c(0,-0.5),depth_test="always",cex=2)
    
  g.xyz.df<-spectrum.df %>% select(X,Y,Z)
  xyzvolume<-convhulln(g.xyz.df)
  triangles3d(g.xyz.df[t(xyzvolume), ], color = rgb(1,1,1), depth_test="always",alpha = 0.4)
  axes3d( edges=c("x--", "y+-", "z--"), col="white", nticks=5,cex=2)
	
})

rgl.snapshot("XYZ.png")




#run this for the CIE xyY chromaticity diagram with embdedded gamuts

g.xy.df<-spectrum.df %>% select(x,y,Y,R,G,B,wl.labels)
g.xy.df

sRGB.df<-data.frame(x=c(0.64,0.3000,0.1500),y=c(0.3300,0.6000,0.0600),R=c(1,0,0),G=c(0,1,0),B=c(0,0,1))
sRGB.df


#Inks
#http://www.color.org/FOGRA39.txt
#INSTRUMENTATION "D50, 2 degree, geometry 45/0, no polarisation filter, white backing, according to ISO 13655"
#C       M     Y       K     X      Y       Z
#100     0     0       0    15.02   22.93   52.85
#100     100   0       0     5.67    4.10   15.67
#0       100   0       0    33.03   16.79   15.01
#0       100   100     0    30.20   16.02    2.30
#0       0     100     0    69.17   74.16    7.04
#100     0     100     0     8.16   18.42    6.74


xy2XYZ<-function(x,y,Y=1) {
	X<-x*Y/y
	Z<-(1-x-y)*Y/y
	c(X,Y,Z)
}

XYZ2xy<-function(X,Y,Z) {
  denom<-X+Y+Z
	x<-X/denom
	y<-Y/denom
	c(x,y)
}

ink.XYZ.scale<-20

ink.df<-data.frame(X=c(15.02,5.67,33.03,30.20,69.17,8.16)/ink.XYZ.scale,
                   Y=c(22.93,4.10,16.79,16.02,74.16,18.42)/ink.XYZ.scale,
                   Z=c(52.85,15.67,15.01,2.30,7.04,6.74)/ink.XYZ.scale)

inks<-as.list(data.frame(apply(ink.df,1,c)))                  
ink.rgb.df<-data.frame(do.call('rbind',lapply(inks,function(x){pmin(1,pmax(0,xyzrgbt %*% x))})))                  
colnames(ink.rgb.df)<-c('R','G','B')
ink.xy.df<-data.frame(do.call("rbind",lapply(inks, function(x){do.call("XYZ2xy",as.list(x))})))
colnames(ink.xy.df)<-c('x','y')
ink.xyz.df<-(t(apply(ink.df,1,function(x){x/sum(x)})))
colnames(ink.xyz.df)=c('x','y','z')
ink.df<-cbind(ink.df,ink.xy.df,ink.rgb.df)
ink.df


#Pigments
#LAB to XYZ: http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_Lab.html

#https://en.wikipedia.org/wiki/Illuminant_D65
xy.refwhite<-c(0.31271,0.32902)
XYZ.refwhite<-c(95.047,100.00,108.883)/25
XYZ.refwhite.list<-as.list(XYZ.refwhite)


#https://people.csail.mit.edu/jaffer/Color/winsor-newton-lab.txt
#for white https://www.gamblincolors.com/wp-content/uploads/2016/01/L-chart-of-all-27-samples-light-aged.pdf
#ivory black dark		CIELab:21.87/0.73/5.13
pigments<-list(yellowochre.Lab=c(67.01,15.24,76.01), cadmiumred.Lab=c(33.1,43.43,36.29), ultramarineblue.Lab=c(30.85,17.31,-72.01))

Lab2XYZ<-function(L,a,b,Xw,Yw,Zw){
	fy<-(L+16)/116
	fz<-fy-(b/200)
	fx<-(a/500)+fy
	
	cie.eps<-216/24389
	cie.kappa<-24389/27
	
	xr<-fx^3
	if(fx^3<=cie.eps){
		xr<-116*(fx^3-16)/cie.kappa
	}
	
	yr<-((L+16)/116)^3
	if(L<=(cie.kappa*cie.eps)){		
		yr<-L/cie.kappa
	}
	
	zr<-fz^3
	if(fz<cie.eps){
		zr<-((116*fz-16))/cie.kappa
	}
	
	c(xr,yr,zr)*c(Xw,Yw,Zw)
}

pigment.xyz.df<-data.frame(do.call("rbind",lapply(pigments, function(x){
	do.call("Lab2XYZ",append(x,XYZ.refwhite.list))
	})))
colnames(pigment.xyz.df)<-c('X','Y','Z')
pigment.xyz.df

pigment.rgb.df<-data.frame(do.call('rbind',lapply(as.list(data.frame(apply(pigment.xyz.df,1,c))) ,function(x){pmin(1,pmax(0,xyzrgbt %*% x))})))                  
colnames(pigment.rgb.df)<-c('R','G','B')
pigment.rgb.df

pigment.xy.df<-data.frame(t(apply(pigment.xyz.df,1,function(x){x/sum(x)})))
colnames(pigment.xy.df)<-c('x','y','z')
pigment.xy.df

pigment.df<-cbind(pigment.xyz.df,pigment.xy.df,pigment.rgb.df)
pigment.df<-pigment.df %>% select(x,y,R,G,B)

label.size<-4.7
dot.size<-4.5



#from https://github.com/colour-science/colour/blob/83b2c224a18fa3a993c37e247c8095b77badb54b/colour/models/dataset/pointer_gamut.py
pointer.gamut.xy.df<-as.data.frame(matrix(c(
  0.659, 0.316,
  0.634, 0.351,
  0.594, 0.391,
  0.557, 0.427,
  0.523, 0.462,
  0.482, 0.491,
  0.444, 0.515,
  0.409, 0.546,
  0.371, 0.558,
  0.332, 0.573,
  0.288, 0.584,
  0.242, 0.576,
  0.202, 0.530,
  0.177, 0.454,
  0.151, 0.389,
  0.151, 0.330,
  0.162, 0.295,
  0.157, 0.266,
  0.159, 0.245,
  0.142, 0.214,
  0.141, 0.195,
  0.129, 0.168,
  0.138, 0.141,
  0.145, 0.129,
  0.145, 0.106,
  0.161, 0.094,
  0.188, 0.084,
  0.252, 0.104,
  0.324, 0.127,
  0.393, 0.165,
  0.451, 0.199,
  0.508, 0.226),
  ncol=2,byrow=T
))
colnames(pointer.gamut.xy.df)<-c('x','y')
pointer.gamut.xy.df

max.xyz.matrix<-matrix(c((spectrum.df %>% select(X,Y,Z,wl.labels) %>% filter(X==max(X)))[1,1],0,0,
                     0,(spectrum.df %>% select(X,Y,Z,wl.labels) %>% filter(Y==max(Y)))[1,2],0,
                     0,0,(spectrum.df %>% select(X,Y,Z,wl.labels) %>% filter(Z==max(Z)))[1,3]), ncol=3, byrow=T)

pri.xyz.xy.df<-as.data.frame(t(apply(max.xyz.matrix,
                       MARGIN=1,
                       FUN=function(x){do.call(XYZ2xy,as.list(t(x)))})))

colnames(pri.xyz.xy.df)<-c("x","y")
pri.xyz.xy.df

pri.lms.xy.df<-as.data.frame(apply( t(lmst),
                                     MARGIN=1,
                                     FUN=function(x){do.call(XYZ2xy,as.list(t(x)))}))


colnames(pri.lms.xy.df)<-c('x','y')
pri.lms.xy.df

ggplot(g.xy.df,aes(x,y,fill="white"))+
geom_polygon(color="grey80",fill="grey80",alpha=0.9)+
scale_color_identity()+
scale_fill_identity()+
scale_shape_identity()+
#geom_polygon(data=pri.xyz.xy.df, aes(x,y),color="white",alpha=0)+
#geom_polygon(data=pri.lms.xy.df, aes(x,y),color="white",alpha=0)+
#geom_polygon(data=pointer.gamut.xy.df, aes(x,y),color="black",fill="grey60",alpha=0.8)+
geom_polygon(data=sRGB.df,color="black",fill="grey40",alpha=0.8)+
geom_point(data=sRGB.df,aes(x,y,color=rgb(R,G,B)), size=dot.size, shape=15)+
geom_polygon(data=ink.df,color="black",fill="grey10",alpha=0.8)+
geom_point(data=ink.df,aes(x,y,color=rgb(R,G,B)), size=dot.size, shape=17)+
#geom_polygon(data=pigment.df,color="black",linetype = 2,fill="grey50",alpha=0.8)+
geom_point(data=pigment.df,aes(x,y,color=rgb(R,G,B)), size=dot.size, shape=3)+
geom_text(data=g.xy.df %>% filter(wl.labels >= 525), aes(x,y,color=rgb(R,G,B),label=wl.labels),hjust=-0.5, vjust=0, cex=label.size)+
geom_text(data=g.xy.df %>% filter(!(wl.labels >= 525), wl.labels>450), aes(x,y,color=rgb(R,G,B),label=wl.labels),hjust=1.3, vjust=0, cex=label.size)+
geom_point(aes(color=rgb(R,G,B)),size=dot.size, shape=16)+
theme(
    axis.text = element_text(color="black",size = 12),
    axis.title=element_text(face="bold", size=14),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "black")
  )
  

ggsave("xyY.png")


