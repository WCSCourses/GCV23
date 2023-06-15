# CUSTOM FUNCTIONS FOR MOVIE PLOTS - uses BEAST mcc trees - you must use read_MCC_tree.R
# S J Lycett
# this version 30 Jun 2022


library(maps)
library(mapdata)
library(mapproj)

# needed to fit HPDs to ellipses
library(conicfit)

# 26 Nov 2021 - added some na.rm=TRUE
fit_HPDs_to_standard <- function( tr=tr, npts=50, ltol=0.005 ) {
  # fit hpds to standard
  hpds <- vector("list",length(tr$lat80))
  for (j in 1:length(tr$lat80)) {
    #xc    <- mean(tr$lat80[[j]])
    #yc    <- mean(tr$lon80[[j]])
    ytemp  <- tr$lat80[[j]] #- xc #tr$latlon[j,1]
    xtemp  <- tr$lon80[[j]] #- yc #tr$latlon[j,2]
    
    xgood  <- (length(xtemp)>3)
    ygood  <- (length(ytemp)>3)
    if (xgood) xgood <- xgood & ((max(xtemp,na.rm=TRUE)-min(xtemp,na.rm=TRUE))>ltol)
    if (ygood) ygood <- ygood & ((max(ytemp,na.rm=TRUE)-min(ytemp,na.rm=TRUE))>ltol)
    
    if (xgood & ygood) {
      kk     <- which(is.finite(xtemp) & is.finite(ytemp))
      xy		 <- cbind(xtemp[kk],ytemp[kk])
      ellipDirect  <- EllipseDirectFit(xy)
      ellipDirectG <- AtoG(ellipDirect)$ParG
    } else {
      if (!xgood) {
        if (!ygood) {
          # small circle only
          ellipDirectG<- c(mean(xtemp,na.rm=TRUE),mean(ytemp,na.rm=TRUE),ltol,ltol,0)
        } else {
          # has y range only
          yrange <- (max(ytemp,na.rm=TRUE)-min(ytemp,na.rm=TRUE))/2
          ellipDirectG<- c(mean(xtemp,na.rm=TRUE),mean(ytemp,na.rm=TRUE),ltol,yrange,0)
        }
      } else {
        # has x range only
        xrange <- (max(xtemp,na.rm=TRUE)-min(xtemp,na.rm=TRUE))/2
        ellipDirectG<- c(mean(xtemp,na.rm=TRUE),mean(ytemp,na.rm=TRUE),xrange,ltol,0)
      }
    }
    
    
    xyFitted    <- calculateEllipse(ellipDirectG[1], ellipDirectG[2], 
                                    ellipDirectG[3], ellipDirectG[4], 
                                    180/pi*ellipDirectG[5], steps=npts)
    
    
    #plot(xyFitted[,1],xyFitted[,2], type="l", col="magenta")
    #points(mean(xtemp),mean(ytemp),pch=21,bg="red")
    #lines( c(mean(xtemp),mean(xtemp)), c(min(ytemp),max(ytemp)), col="blue")
    #lines( c(min(xtemp),max(xtemp)), c(mean(ytemp),mean(ytemp)), col="blue")
    #if (length(xtemp)==length(ytemp) & length(xtemp)>1) {
    #	polygon(xtemp,ytemp)			
    #} else {
    #	if (length(xtemp)==1) {
    #		points( array(xtemp, length(ytemp)), ytemp)
    #	}
    #	if (length(ytemp)==1) {
    #		points( xtemp, array(ytemp, length(xtemp)) )
    #	}
    #}
    
    hpds[[j]]	<- xyFitted
  }
  tr$hpds <- hpds
  
  # hpds differences
  hpd_diffs <- vector("list", length(tr$edge.length))
  for (j in 1:length(tr$edge.length)) {
    fromHPD <- hpds[[ tr$edge[j,1] ]]
    toHPD   <- hpds[[ tr$edge[j,2] ]]
    hpd_diffs[[j]] <- toHPD-fromHPD
  }
  tr$hpd_diffs <- hpd_diffs
  
  
  return( tr )
}


# 17 mar 2019
addFullTraitSet <- function(tr=tr, propIndex=1) {
  # see read_MCC_tree.R for getFullTraitSet function
  fullset 	<- getFullTraitSet(tr=tr, propIndex=propIndex)
  tr$fullset 	<- fullset
  return( tr )
}

# 17 mar 2019
interpolate_discrete_element <- function( fromSet=fromSet, toSet=toSet, fractTime=fractTime) {
  ff 	 <- matrix( rep(fractTime,length(fromSet[1,])), length(fractTime), length(fromSet[1,]) )
  midSet <- (fromSet*(1-ff)) + (toSet*ff)
  return( midSet )
}

# returns nodes and hpds suitable for plotting at a time point
# interpolates between branches which cut the time point
# optionally return the interpolated discrete trait if addFullTraitSet is used
pts_and_hpds_at_time <- function(timePt, tr=tr, ii=c(), xlim=c(-180,180), ylim=c(-60,80)) {
  
  # latlon -2 = x, -1 = y
  fromY		<- tr$latlon[tr$edge[,1],1]
  fromX 	<- tr$latlon[tr$edge[,1],2]
  toY		<- tr$latlon[tr$edge[,2],1]
  toX 		<- tr$latlon[tr$edge[,2],2]
  fromTime 	<- tr$nodeTimes[tr$edge[,1]]
  toTime   	<- tr$nodeTimes[tr$edge[,2]]
  branchTime  <- (toTime-fromTime)
  branchX	<- (toX-fromX)
  branchY	<- (toY-fromY)
  
  # 16 mar 2019
  hpds 		<- tr$hpds
  hpd_diffs 	<- tr$hpd_diffs
  
  # identify branches
  if (length(ii)==0) {
    ii 		<- which( toTime >= timePt & fromTime <= timePt  &
                     toX >= xlim[1] & toX <= xlim[2] & toY >= ylim[1] & toY <= ylim[2] )
  }
  
  fractTime 	<- 1-(toTime-timePt)/(branchTime)
  
  if (length(ii)>0) {
    if (length(ii)==1) {
      ii <- c(ii,ii)
    }
    # get mcc pts
    x_pos		<- (fractTime*branchX + fromX)[ii]
    y_pos		<- (fractTime*branchY + fromY)[ii]
    
    pts 		<- list(x=x_pos,y=y_pos)
    
    # hpds
    pts_polys <- vector("list",length(ii))
    for (j in 1:length(ii)) {
      fromNode	  <- tr$edge[ii[j],1]
      temp_hpd_diff <- hpd_diffs[[ii[j]]]
      temp_hpd	  <- fractTime[ii[j]]*temp_hpd_diff + hpds[[fromNode]]
      pts_polys[[j]] <- temp_hpd
    }
  } else {
    pts       <- NULL
    pts_polys <- NULL
  }
  
  if ( any(attributes(tr)$names=="fullset") ) {
    # interpolate the selected discrete trait
    # use for colouring
    fromSet <- tr$fullset[tr$edge[ii,1],]
    toSet	  <- tr$fullset[tr$edge[ii,2],]
    midSet  <- interpolate_discrete_element(fromSet=fromSet, toSet=toSet, fractTime=fractTime[ii])
    interpol_trait <- colnames(midSet)[apply(midSet, 1, which.max)]
    
    return( list(pts=pts, pts_polys=pts_polys, 
                 ii=ii, fractTime=fractTime[ii], 
                 interpol_trait=interpol_trait, midSet=midSet) )
  } else {
    return( list(pts=pts, pts_polys=pts_polys, ii=ii, fractTime=fractTime[ii]) )
  }
}

# 17 mar 2019; 5 Nov 2021
plot_at_time <- function(tpt=tpt, timePt=timePt, tr=tr, xlim=c(-180,180), ylim=c(-60,80), 
                         show.hpds=TRUE, solid.pts=TRUE, show.legend=TRUE, new.plot=TRUE,
                         timelegpos="bottomright",timeFormat="ddmmyy",
                         legpos="bottomleft",
                         fcol="white", bgcol="grey90", bdcol="grey70", fill=TRUE,
                         use.fullset.cols=TRUE,
                         h1=0, h2=0,
                         s1=0.3, s2=0.7,
                         b1=0.9, b2=0.7,
                         t1=0.25, t2=0.75, useWorldHires=TRUE ) {
  
  if ( length(tpt)== 0 ) {
    tpt 	<- pts_and_hpds_at_time(timePt, tr=tr, xlim=xlim, ylim=ylim)
  }
  
  pcol 	<- array(hsv(h1, s1, b1, t1), length(tpt$pts_polys))
  pcol2 <- array(hsv(h2, s2, b2, t2), length(tpt$pts_polys))
  
  if (new.plot) {
    if (useWorldHires) {
      map("worldHires", xlim=xlim, ylim=ylim, 
          col=fcol, fill=fill, border=bdcol, bg=bgcol)
    } else {
      map("world", xlim=xlim, ylim=ylim, 
          col=fcol, fill=fill, border=bdcol, bg=bgcol)
    }
  }
  
  if ( use.fullset.cols & any(attributes(tr)$names=="fullset") ) {
    utraits <- colnames(tr$fullset)
    tcols   <- get_BEAST_cols( length(utraits), sat=s1, bright=b1, transparency=t1 )
    tcols2  <- get_BEAST_cols( length(utraits), sat=s2, bright=b2, transparency=t2 )
    
    tinds	  <- match(tpt$interpol_trait, utraits)
    pcol	  <- tcols[tinds]
    pcol2	  <- tcols2[tinds]
    
    if (show.legend) {
      legend(legpos, paste(utraits), pch=21, pt.bg=tcols2, bty="n")
    }
  }
  
  if (show.hpds) {
    for (j in 1:length(tpt$pts_polys)) {
      polygon(tpt$pts_polys[[j]], col=pcol[j], bg=pcol[j], border=FALSE)
    }
  }
  if (solid.pts) {
    points(tpt$pts, pch=21, bg=pcol2)
  } else {
    points(tpt$pts, col=pcol2)
  }
  
  if (show.legend) {
    if (timeFormat=="ddmmyy" | timeFormat=="mmyy") {
      timeTxt <- invertDecimalDate(timePt,ddmmyy=TRUE)
      if (timeFormat=="mmyy") {
        timeTxt <- substring(timeTxt, 4)
      }
      legend(timelegpos, timeTxt, bty="n", pch=NA)
    } else {
      legend(timelegpos, format(timePt,digits=7), bty="n", pch=NA)
    }
  }
}


plot_mcc_tree_with_hpds <- function(tr=tr, xlim=c(-180,180), ylim=c(-60,80), 
                                    tlim=c(NA,NA),
                                    propIndex=0,
                                    show.hpds=TRUE, solid.pts=TRUE, show.legend=TRUE, 
                                    timelegpos="bottomright",
                                    legpos="bottomleft",
                                    fcol="white", bgcol="grey90", bdcol="grey70", fill=TRUE,
                                    s1=0.3, s2=0.7,
                                    b1=0.9, b2=0.7,
                                    t1=0.25, t2=0.75) {
  
  pcol 	<- array(hsv(0, s1, b1, t1), length(tr$hpds))
  pcol2 <- array(hsv(0.6, s2, b2, t2), length(tr$hpds))
  ecol2 <- array(hsv(0, 0, 0.2, t2), length(tr$edge))
  
  map("worldHires", xlim=xlim, ylim=ylim, 
      col=fcol, fill=fill, border=bdcol, bg=bgcol)
  
  if (propIndex>0) {
    utraits <- tr$uprops[[propIndex]]
    tcols   <- get_BEAST_cols( length(utraits), sat=s1, bright=b1, transparency=t1 )
    tcols2  <- get_BEAST_cols( length(utraits), sat=s2, bright=b2, transparency=t2 )
    
    tinds	  <- match(tr$props[,propIndex], utraits)
    pcol	  <- tcols[tinds]
    pcol2	  <- tcols2[tinds]
    einds	  <- match(tr$props[tr$edge[,1],propIndex], utraits)
    ecol2	  <- tcols2[einds]
    
    if (show.legend) {
      legend(legpos, paste(utraits), pch=21, pt.bg=tcols2, bty="n")
    }
  }
  
  tr$ntips <- length(tr$tip.label)
  tinds    <- 1:tr$ntips
  ninds    <- (tr$ntips+1):length(tr$latlon[,1])
  fromY		<- tr$latlon[tr$edge[,1],1]
  fromX 	<- tr$latlon[tr$edge[,1],2]
  toY		<- tr$latlon[tr$edge[,2],1]
  toX 		<- tr$latlon[tr$edge[,2],2]
  
  if (!is.na(tlim[1])) {
    ok_nodes <- which(tr$nodeTimes >= tlim[1] & tr$nodeTimes <= tlim[2])
    tinds	   <- intersect(tinds, ok_nodes)
    ninds	   <- intersect(ninds, ok_nodes)
    ok_edges <- which(tr$nodeTimes[tr$edge[,2]] >= tlim[1] & 
                        tr$nodeTimes[tr$edge[,1]] <= tlim[2] )
  } else {
    ok_nodes <- c(tinds,ninds)
    ok_edges <- 1:length(tr$edge[,1])
  }
  
  #for (i in 1:length(tr$hpds)) {
  for (j in length(ok_nodes):1) {
    i <- ok_nodes[j]
    polygon(tr$hpds[[i]], col=pcol[i], bg=pcol[i], border=FALSE)
  }
  
  arrows( 	fromX[ok_edges], 	fromY[ok_edges], 
           toX[ok_edges], 	toY[ok_edges], 
           length=0.1, angle=15, col=ecol2[ok_edges])
  
  points(tr$latlon[tinds,2], tr$latlon[tinds,1], pch=21, bg=pcol2[tinds])
  points(tr$latlon[ninds,2], tr$latlon[ninds,1], pch=21, col=pcol2[ninds])
  
}
