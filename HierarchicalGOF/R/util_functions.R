#' plot a map of estimated abundance or related quantity on a grid.  
#' @param N quantity to plot
#' @param Coords data.frame holding x and y coordinates of each prediction cell
#' @param highlight [optional] A vector of cells to highlight.  Default = NULL
#' @param myPalette A prespecified palette from RColorBrewer (default is YlOrBr(9) )
#' @param leg.title Title of the legend
#' @return ggplot2 map of predictions on a grid
#' @export 
#' @author Paul Conn \email{paul.conn@@noaa.gov}
plot.prediction.map<-function(N,Coords,highlight=NULL,myPalette=NULL,leg.title="Abundance"){
  #require(rgeos)
  #require(ggplot2)
  #library(RColorBrewer)
  if(is.null(myPalette))myPalette <- colorRampPalette(brewer.pal(9, "YlOrBr"))
  if(is.null(highlight)==FALSE){
    midpoints=Coords[highlight,]
    colnames(midpoints)=c("Easting","Northing")
  }
  Abundance=N
  Cur.df=cbind(Coords,Abundance)
  new.colnames=colnames(Cur.df)
  new.colnames[1:2]=c("Easting","Northing")
  colnames(Cur.df)=new.colnames
  Grid.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
  p1=ggplot(Cur.df)+aes(Easting,Northing,fill=Abundance)+geom_raster()+Grid.theme+scale_fill_gradientn(colours=myPalette(100),name=leg.title)
  if(is.null(highlight)==FALSE){
    #p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067,xmax=Easting,ymin=Northing,ymax=Northing+25067))
    p1=p1+geom_rect(data=midpoints,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
  }
  p1
}

