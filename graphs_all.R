library(ggplot2)
library(dplyr)
library(tidyr)
library(colorscience)
library(ggspectra)
library(scatterplot3d)
library(rgl)
library(geometry)

o.df<-StockmanMacLeodJohnson2degCIEadjConeFundamentals1993

colnames(o.df)=c("wavelength","L","M","S")
o.df[,c('L','M','S')]<-10^(o.df[,c('L','M','S')])
o.df

purple1.wl<-420
purple2.wl<-630
lms.df <- o.df %>% filter(wavelength>=purple1.wl,wavelength<=purple2.wl )
lms.df

lms.purple1<-lms.df %>% filter(wavelength==purple1.wl)
lms.purple2<-lms.df %>% filter(wavelength==purple2.wl)

purple.funcs<-apply(rbind(lms.purple1,lms.purple2), 2, approxfun, x=c(purple1.wl, purple2.wl))
purples.lms.df<-data.frame(t(sapply(seq(purple1.wl,purple2.wl,25),function(x){c(NA,purple.funcs$L(x),purple.funcs$M(x),purple.funcs$S(x))})))
colnames(purples.lms.df)=colnames(lms.df)

lms.df<-rbind(lms.df,purples.lms.df)

lmss<-as.list(data.frame(apply(lms.df[,c("L","M","S")],1,c)))
names(lmss)<-lms.df$wavelength
lmss

#von Kries from wikipedia https://en.wikipedia.org/wiki/LMS_color_space
#lmst<-solve(t(matrix(c( 0.38971, 0.68898, -0.07868,
#                       -0.22981, 1.18340, 0.04641,
#                        0,       0,        1), 3,3)))

#D65 norm'd
lmst<-solve(t(matrix(c( 0.4002, 0.7076,   -0.0808,
                       -0.2263, 1.1653,    0.0457,
                        0,       0,        0.9182), 3,3)))

                              

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


xyzrgbt

lms2xyz<-function(x){lmst %*% x}
                              
lms2rgbtt<-function(x){((xyzrgbt %*% lmst) %*% x)}

xyz.df<-data.frame(do.call("rbind",lapply(lmss,function(x){t(lms2xyz(x))})))
colnames(xyz.df)<-c('X','Y','Z')


#rgb.df<-do.call("rbind",lapply(lmss,function(x){pmin(1,pmax(0,lms2rgbtt(x)))}))
rgb.df<-do.call("rbind",lapply(lmss,function(x){pmin(1,pmax(0,sapply(lms2rgbtt(x),gammacor)))}))
colnames(rgb.df)<-c('R','G','B')

spectrum.df<-cbind(lms.df,xyz.df,rgb.df)
lmsnorm.df<-spectrum.df %>% select(Ln=L,Mn=M,Sn=S)
chromaticity.df<-spectrum.df %>% select(x=X,y=Y,z=Z)

lmsnorm.df<-t(apply(lmsnorm.df,1,function(x){x/sum(x)}))
chromaticity.df<-t(apply(chromaticity.df,1,function(x){x/sum(x)}))

spectrum.df<-cbind(spectrum.df,lmsnorm.df,chromaticity.df) %>% mutate(wl.labels=as.character(wavelength))
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



#cone spectral sensitivities
t.df <- o.df %>% gather(cone.type, response, L:S)

ggplot(t.df, aes(x=wavelength,y=response,group=cone.type))+
stat_wl_strip(alpha = 1, ymin = -Inf, ymax = Inf)+
scale_fill_identity()+
geom_line(color="white",size=2)+
geom_text(data=data.frame(cone.type=c('L','M','S'),wavelength=c(605,505,460),response=c(0.8,0.8,0.8)), aes(label=cone.type),size=8,hjust=0,vjust=-1,color="white")+
labs(x="wavelength (nm)", y="normalized cone sensitivity")+
xlim(400,680)+
theme(
    axis.text = element_text(color="black",size = 12),
    axis.title=element_text(face="bold", size=14),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "black")
  )


ggsave("cones.png")


#run this 'with' block for the 3D LMS chart
with(spectrum.df, {
 	if( rgl.cur() == 0 ) {
    rgl.open()
    par3d(zoom=0.8227028, userMatrix=matrix(c(0.972281754016876,0.199476152658463,0.121972970664501,0,0,0.521670699119568,-0.853146910667419,0,-0.23381219804287,0.829499185085297,0.507210910320282,0,-0.0236389161640957,0.175143327954142,2.32254899756867e-08,1),4,4),windowRect = c( 1440, 45, 3440, 1441 ), antialias=4 )
    rgl.bg(color = "black" )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
 # rgl.viewpoint(theta = 0, phi = -130, zoom = 1)
  axes.title.size=2
  axes3d( edges=c("x--", "y+-", "z--"), col="white", cex=1.5,nticks=5)
	
	box.color=rgb(0.6,0.6,0.6)
	box.alpha=0.5
	outset=rep(-0.01,4)

    rgl.lines(c(0,0),c(0,1),c(0,0),color="white")
    rgl.lines(c(0,0),c(1,1),c(0,1),color="white")
    rgl.lines(c(0,1),c(1,1),c(0,0),color="white")
    	
    
    rgl.texts(rbind(c(0.5,0,0),c(0,0,0.5)), text = c("L", "S"), color = "white",
             adj = c(8, 8), cex = axes.title.size)
    rgl.texts(rbind(c(1,0.5,0)), text = c("M"), color = "white",
             adj = c(-5, 0), cex = axes.title.size)        
    rgl.points(L,M,S,col=rgb(spectrum.df %>% select(R,G,B)),size=6) 
    text3d(L,M,S,wl.labels,adj=c(0,0),depth_test="always", cex=1.5)
    
    g.lms.df<-spectrum.df %>% select(L,M,S)
    lmsvolume<-convhulln(g.lms.df)
    triangles3d(g.lms.df[t(lmsvolume), ], color = rgb(1,1,1), depth_test="always",alpha = 0.4)
})

rgl.snapshot("LMS.png")




#run this 'with' block for the 3D XYZ chart
with(spectrum.df, {
	if( rgl.cur() == 0 ) {
    rgl.open()
   par3d(zoom=0.8227028, userMatrix=matrix(c(0.972281754016876,0.199476152658463,0.121972970664501,0,0,0.521670699119568,-0.853146910667419,0,-0.23381219804287,0.829499185085297,0.507210910320282,0,-0.0236389161640957,0.175143327954142,2.32254899756867e-08,1),4,4),windowRect = c( 1440, 45, 3440, 1441 ), antialias=4)
    rgl.bg(color = rgb(0,0,0) )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
 	
  box.color=rgb(0.6,0.6,0.6)
  box.alpha=0.5
  outset=rep(-0.01,4)

  rgl.lines(c(0,0),c(0,1),c(0,0),color="white")
  rgl.lines(c(0,0),c(1,1),c(0,1),color="white")
  rgl.lines(c(0,1),c(1,1),c(0,0),color="white")
    	
  axes.title.size=2
  rgl.texts(rbind(c(0.5,0,0),c(0,0,0.5)), text = c("X", "Z"), color = "white",
             adj = c(8, 8), cex = axes.title.size)
  rgl.texts(rbind(c(1,0.5,0)), text = c("Y"), color = "white",
             adj = c(-10, 0), cex = axes.title.size)        
  rgl.points(X,Y,Z,col=rgb(spectrum.df %>% select(R,G,B)),size=6,depth_test="always") 
  text3d(X,Y,Z,wl.labels,adj=c(0,0),depth_test="always",cex=1.5)
    
  g.xyz.df<-spectrum.df %>% select(X,Y,Z)
  xyzvolume<-convhulln(g.xyz.df)
  triangles3d(g.xyz.df[t(xyzvolume), ], color = rgb(1,1,1), depth_test="always",alpha = 0.4)
  axes3d( edges=c("x--", "y+-", "z--"), col="white", nticks=5,cex=1.5)
	
})

rgl.snapshot("XYZ.png")




#run this for the CIE xyY chromaticity diagram with embdedded gamuts

g.xy.df<-spectrum.df %>% select(x,y,Y,R,G,B,wavelength,wl.labels)
g.xy.df

sRGB.df<-data.frame(x=c(0.64,0.3000,0.1500),y=c(0.3300,0.6000,0.0600),R=c(1,0,0),G=c(0,1,0),B=c(0,0,1))
sRGB.df


#Inks
#http://www.color.org/FOGRA39.txt
#INSTRUMENTATION "D50, 2 degree, geometry 45/0, no polarisation filter, white backing, according to ISO 13655"
#C       M     Y       K     X      Y       Z
#100     0     0       0    15.02   22.93   52.85
#0       100   0       0    33.03   16.79   15.01
#0       0     100     0    69.17   74.16    7.04

fcmyk.xy.df<-data.frame(x=c(15.02,33.03,69.17),y=c(22.93,16.79,74.16),z=c(52.85,15.01,7.04))
fcmyk.xy.df<-data.frame(t(apply(fcmyk.xy.df,1,function(x){x/sum(x)})))
fcmyk.xy.df<-fcmyk.xy.df %>% mutate (R=c(0,1,1),G=c(1,0,1),B=c(1,1,0))




#Pigments
#LAB to XYZ: http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_Lab.html
#https://www.goldenpaints.com/products/colors/heavy-body/yellow-oxide  L*61.2 a*15.03 b*43.46
#https://www.goldenpaints.com/products/colors/heavy-body/c-p--cadmium-red-dark  L*38.2 a*42.27 b*19.41
#https://www.goldenpaints.com/products/colors/heavy-body/ultramarine-blue  L*24.36 a*12.28 b*-25.74


#https://en.wikipedia.org/wiki/Illuminant_D65
xy.refwhite<-c(0.31271,0.32902)
XYZ.refwhite<-c(95.047,100.00,108.883)/100
XYZ.refwhite.list<-as.list(XYZ.refwhite)

yellowochre.Lab<-c(61.2,15.03,43.46)
cadmiumred.Lab<-c(38.2,42.47,19.41)
ultramarineblue.Lab<-c(24.36,12.28,-25.74)


Lab2XYZ<-function(L,a,b,Xw,Yw,Zw){
	fy<-(L+16)/116
	fz<-fy-(b/200)
	fx<-(a/500)+fy
	
	cie.eps<-216/24389
	cie.kappa<-24389/27
	
	xr<-	fx^3
	if(fx^3<=cie.eps){
		xr<-116*(fx^3-16)/cie.kappa
	}
	
	yr<-((L+16)/116)^3
	if(L<=cie.kappa*cie.eps){		
		yr<-L/cie.kappa
	}
	
	zr<-fz^3
	if(fz<cie.eps){
		zr<-((116*fz-16))/cie.kappa
	}
	
	c(xr,yr,zr)*c(Xw,Yw,Zw)
}

pigment.xy.df<-data.frame(do.call("rbind",lapply(list(ub=ultramarineblue.Lab,cr=cadmiumred.Lab,yo=yellowochre.Lab), function(x){
	do.call("Lab2XYZ",append(x,XYZ.refwhite.list))
	})))

pigment.xy.df

colnames(pigment.xy.df)<-c('x','y','z')
pigment.xy.df<-data.frame(t(apply(pigment.xyz.df,1,function(x){x/sum(x)})))
pigment.xy.df<-pigment.xy.df %>% select(x,y)
pigment.xy.df<-pigment.xy.df %>% mutate (R=c(0,1,1),G=c(0,0,1),B=c(1,0,0))

ggplot(g.xy.df,aes(x,y,fill="white"))+
geom_point(aes(color=rgb(R,G,B)))+
geom_polygon(color="grey90",alpha=0.5)+
scale_color_identity()+
scale_fill_identity()+
geom_point(data=sRGB.df,aes(x,y,color=rgb(R,G,B)), size=8)+
geom_polygon(data=sRGB.df,color="grey70",alpha=0.5)+
geom_point(data=fcmyk.xy.df,aes(x,y,color=rgb(R,G,B)), size=8)+
geom_polygon(data=fcmyk.xy.df,color="grey50",alpha=0.5)+
geom_point(data=pigment.xy.df,aes(x,y,color=rgb(R,G,B)), size=8)+
geom_polygon(data= pigment.xy.df,color="grey30",alpha=0.5)+
geom_text(data=g.xy.df %>% filter(wavelength >= 525), aes(x,y,color=rgb(R,G,B),label=wl.labels),hjust=-0.5, vjust=0, cex=5)+
geom_text(data=g.xy.df %>% filter(!(wavelength >= 525)), aes(x,y,color=rgb(R,G,B),label=wl.labels),hjust=1, vjust=0, cex=5)+
theme(
    axis.text = element_text(color="black",size = 12),
    axis.title=element_text(face="bold", size=14),
    panel.grid.major = element_line(colour = "grey50"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "black")
  )

ggsave("xyY.png")


