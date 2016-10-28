###############
## Modified from Carter T Butts' original genpop functon
## Builds a list of household population
###############
#hh<-c("hh.1person","hh.2person","hh.3person", "hh.4person","hh.5person","hh.6person")
#hh.size<-6
#total.pop<-"pop2000"
#test<-genpop(minnesota.tract,hh=c("hh.1person","hh.2person","hh.3person", "hh.4person","hh.5person","hh.6person"),total.pop="pop2000",hh.size=6)

genpop<-function(spa,hh,total.pop,hh.size,verbose=FALSE){
	pop<-list()
		for(i in 1:length(spa@polygons)){
			hpop<-sapply(hh,function(x){spa@data[i,x]})
			 bpop<-spa@data[i,total.pop]
			if((hpop%*%(1:hh.size))<bpop){
				if(verbose){
				cat("FYI, down by",bpop-hpop%*%(1:hh.size),"people.  Correcting.\n")
				}
				hpop[1]<-hpop[1]+bpop-hpop%*%(1:hh.size)
				}
		pop[[i]]<-cbind(1:hh.size,hpop)
		}
	pop
}

#genpopIndv(minnesota.tract,"pop2000")
genpopIndv<-function(spa,total.pop){lapply(spa@data[,total.pop],function(x){cbind(1,x)})}




#################
### test<-rspop.SpatialPolygonDataFrame(geo=poly,method="uniform",list.pop=list(hh=c("hh.1person","hh.2person","hh.3person", "hh.4person","hh.5person","hh.6person"),total.pop="pop2000",hh.size=6))
## test2<-rspop.SpatialPolygonDataFrame(geo=poly,method="uniform",list.pop=list(total.pop="pop2000"),household=FALSE)
#################
## Register S3
#S3method(rspop, SpatialPolygonDataFrame)

rspop.SpatialPolygonDataFrame<-function(geo,method = c("uniform", "halton"), stack.rad = 10,stack.dis = 4, household.jitter = 5,longLat=TRUE,household=TRUE,list.pop){
#require(rgdal)
if(household){
	pop<-do.call("genpop",c(list(geo),list.pop))
}else{
	pop<-do.call("genpopIndv",c(list(geo),list.pop))
}

projection.point=coordinates(buildBB(geo))
coord<-rspop(geo,pop,method=method,stack.rad=stack.rad,stack.dis=stack.dis,projection.point=projection.point)

if(longLat){
coordLL<-sp::SpatialPoints(cbind(coord$x,coord$y),proj4string=sp::CRS(paste("+proj=ortho +lat_0=",projection.point[2], " +lon_0=", projection.point[1],collapse = "", sep = "")))
coordLL<-sp::coordinates(sp::spTransform(coordLL,sp::CRS(proj4string(geo))))
coord$x<-coordLL[,1]
coord$y<-coordLL[,2]
}
coord
}





################
## Tools for helping with plotting
################

buildBB<-function(polygon,bb.epsilon=0){
  	temp1 <- sp::bbox(polygon)
        temp1[1, 1] <- temp1[1, 1] - bb.epsilon
        temp1[1, 2] <- temp1[1, 2] + bb.epsilon
        temp1[2, 1] <- temp1[2, 1] - bb.epsilon
        temp1[2, 2] <- temp1[2, 2] + bb.epsilon
        temp <- matrix(c(temp1[1, 1], temp1[2, 1], temp1[1, 2], 
            temp1[2, 1], temp1[1, 2], temp1[2, 2], temp1[1, 1], 
            temp1[2, 2]), ncol = 2, nrow = 4, byrow = TRUE)
        pol<-sp::Polygon(rbind(temp, c(temp[1,1],temp[1,2])), hole = FALSE)
        pol<-sp::Polygons(list(pol),"ID")
        sp::SpatialPolygons(list(pol), proj4string = sp::CRS("+proj=longlat +datum=NAD83"))
  }

orthoProj<-function(polygon,pp){
	#require(rgdal)
	projection.point<-pp
	trans<-sp::CRS(paste("+proj=ortho +lat_0=",projection.point[2], " +lon_0=", projection.point[1],collapse = "", sep = ""))
	rgdal::spTransform(polygon,trans)
}

gplot.spatial<-function(x,coord,bb,list.edges,list.vertex,background=NULL,list.background,list.axis=list(one=list(side=1,cex.axis=.4,lwd=.5,las=1),two=list(side=2,cex.axis=.4,lwd=.5)),axis=TRUE,bg=NULL){
	plot.new()
	plot.window(xlim=bb[1,],ylim=bb[2,],asp=1)

if(!is.null(bg)){
	plot(background,border=bg,col=bg,add=TRUE)
}



	if(axis){
		box()
		do.call("axis",list.axis[[1]])
		do.call("axis",list.axis[[2]])
	}
	
	if(!is.null(background)){do.call("plot",c(list(x=background),list.background,list(add=TRUE)))}
	x<-as.matrix.network(x,"edgelist")
	do.call("lines",c(list(x=rbind(coord[x[,1],],coord[x[,2],])),list.edges))           
	do.call("points",c(list(x=coord),list.vertex))        
}



