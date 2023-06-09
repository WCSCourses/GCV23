# function to calculate decimal date
# S. J. Lycett
# 27 Oct 2010
# 15 Dec 2010 - defaultDay & defaultMonth set to 1, if defaultDay & DefaultMonth = -1 then random day and month are chosen
# from matrixTranslate.R version 3 June 2010

# 19 May 2011
# from matrixTranslate.R version 7 April 2011

# 4 August 2011 - corrected calcDecimalDate_from_yymmdd (defaultMonth was wrong)
# 15 Mar 2013 - calcDecimalDate_from_yymmdd will now run if e.g. 1934//

# 8 May 2020 (actually 7 May 2020 at 5pm) - now with leap years
# 24 June 2020 - still trying to solve rounding errors
# 3 Aug 2020 - still trying to solve round errors
# 8 Aug 2020 - more rounding errors (invertDecimalDate)

calcDecimalDate	<- function(day, month, year, defaultMonth=6, defaultDay=15) {
	# 8 May 2020
	daysInMonth <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
	daysInYear  <- 365
	if ( (year %% 4)== 0 ) {
		daysInMonth[2] <- 29
		daysInYear <- 366
	}
	cd <- array(0,12)
	for (j in 2:12) {
		cd[j] <- cd[j-1]+daysInMonth[j-1]
	}

	#original
	#cd	<- c(0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334)

	if (month==0) {
		if (defaultMonth >= 1) {
			month <- defaultMonth
		} else {
			month	<- ceiling(runif(1)*12)
		}
	}

	if (day==0) {
		if (defaultDay >= 1) {
			day	<- defaultDay
		} else {
			day	<- ceiling(runif(1)*30)
		}
	}

	dd	<- cd[month] + day - 1
	
	decDate <- year + (dd/daysInYear)

	return ( decDate )
}


invertDecimalDate <- function( decDate, formatAsTxt=FALSE, ddmmyy=FALSE ) {
	#cd	<- c(0,  31,  59,  90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
	#fractD<- cd/365
  
  year		<- floor(decDate)
  
  # 8 May 2020
  daysInMonth <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  daysInYear <- 365
  if ( (year %% 4)== 0 ) {
    daysInMonth[2] <- 29
    daysInYear <- 366
  }
  cd <- array(0,12)
  for (j in 2:12) {
    cd[j] <- cd[j-1]+daysInMonth[j-1]
  }
  cd <- c(cd,daysInYear)
  fractD <- cd/daysInYear

	fractYear 	<- decDate-year
	#month		<- which(fractD >= fractYear)[1]-1
	
	# rounding errors causing issues
	# 24 june 2020 - rounding erros still casusing issues
	#tol <- 0.1/daysInYear
	#month <- which( (fractYear-fractD) >= -tol )
	#month <- month[length(month)]
	
	# 24 june 2020- do by days
	tol   <- 0.05
	month <- which(((fractYear*daysInYear) - cd) >= -tol)
	month <- month[length(month)]

	if (month > 0) {
		fractMonth  <- fractYear-fractD[month]
		day		<- round((fractMonth*daysInYear)+1)
		
		#8 Aug 2020 - correction for just before year boundary (e.g. caused by 2020 - 1e-9..)
		if (month==13) {
		  year <- year+1
		  month<- 1
		}
		
		if (day > daysInMonth[month]) {
		  day   <- 1
		  month <- month+1
		  if (month>12) {
		    month <- 1
		    year  <- year+1
		  }
		}
	} else {
	  month <- 1
		day   <- 1
	}

	res <- c(day,month,year)

	if (formatAsTxt) {
		if (month < 10) {
			mm  <- paste("0",month,sep="")
		} else {
			mm <- month
		}
	  if (day < 10) {
	    dd <- paste("0",day,sep="")
	  } else {
	    dd <- paste(day)
	  }
		res <- paste(year,mm,dd,sep="-")
	}

	if (ddmmyy) {
		months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
		if (day < 10) {
			dd <- paste("0",day,sep="")
		} else {
			dd <- day
		}
		res <- paste(dd,months[month],year,sep="-")
	}
	return( res )
	
}


# 4 Nov 2013 - change defaults
# 10 Dec 2014 - if only year then -> year + 0.5, else 15th day of month (Good for NCBI flu)
# 5 Feb 2016 - named months
calcDecimalDate_fromTxt	<- function( dateTxt, sep="/", namedMonths=FALSE, dayFirst=FALSE) {
	els 	<- strsplit(dateTxt, sep)[[1]]
	if (dayFirst) {
		if (length(els) > 1) {
			els <- els[length(els):1]
		}
	}

	year 	<- as.integer(els[1])

	if (length(els)==1) {
		month <- 6  #7
		day	<- 15 #2
		decDate <- year + 0.5
	} else {
	
		if (length(els)==2) {
			if (nchar(els[2]) > 0) {
				if (namedMonths) {
					month <- match(els[2], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
				} else {
					month <- as.integer(els[2])
				}
				day	<- 15
				decDate <- calcDecimalDate(day,month,year)
			} else {
				month <- 6 #7
				day   <- 15 #2
				decDate <- year + 0.5
			}
		} else {
			if (namedMonths) {
				month <- match(els[2], c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
			} else {
				month <- as.integer(els[2])
			}

			if (nchar(els[3]) > 0) {
				day 	<- as.integer(els[3])
			} else {
				day <- 15
			}
			decDate <- calcDecimalDate(day,month,year)
		}
	}


	return ( decDate )
}


calcDecimalDate_from_yymmdd	<- function( dateTxt, sep="/", ycutoff=15, defaultMonth=6, defaultDay=15 ) {
	els	<- strsplit(dateTxt, sep)[[1]]
	yy	<- as.integer(els[1])
	mm	<- as.integer(els[2])
	dd	<- as.integer(els[3])

	if (!is.finite(yy)) {
		return( -1 )
	} else {
		if (yy <= ycutoff) {
			yy <- yy+2000
		}
		if ((yy > ycutoff) & (yy < 99)) {
			yy <- yy+1900
		}

		if (!is.finite(mm)) {
			mm <- 0
		}
		if (!is.finite(dd)) {
			dd <- 0
		}
		return ( calcDecimalDate( dd, mm, yy, defaultMonth=defaultMonth, defaultDay=defaultDay ) )
	}
	
}

########################################################################
# 8 May 2020 - check the decimal dates

doCheck <- FALSE
if (doCheck) {
  
  
  daysInMonth <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  years <- c(2019,2020)
  days  <- c(1,2,29,29,30,31)
  
  first <- TRUE
  
  for (i in 1:length(years)) {
    for (k in 1:length(days)) {
      
      yy <- array(years[i],12)
      mm <- 1:12
      dd <- array(days[k],12)
      for (j in 1:12) {
        
        if (dd[j] > daysInMonth[mm[j]]) {
          ok <- FALSE
        } else {
          ok <- TRUE
        }
        
        decDate <- calcDecimalDate(dd[j],mm[j],yy[j])
        invDate <- invertDecimalDate(decDate)
        
        dd1 <- paste(c(dd[j],mm[j],yy[j]),collapse="-")
        dd2 <- paste(invDate,collapse="-")
        dateMatch <- dd1==dd2
        if (ok) {
          dateOK <- dateMatch
        } else {
          dateOk <- !dateMatch
        }
        temp <- c(dd1,dd2,dateMatch,dateOK)
        if (first) {
          res <- temp
          first<- FALSE
        } else {
          res  <- rbind(res, temp)
        }
      }
      
    }
  }
  
  colnames(res) <- c("OriginalDate","Reinvert","DatesMatch","DatesOK")
  rownames(res) <- NULL
  print(res)
  

}