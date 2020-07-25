#LOWESS Smoothing
#aka LOESS smoothing
#Date 7/18/2020 Jae Chang
#explained in: https://www.itl.nist.gov/div898/handbook/pmd/section1/pmd144.htm

#~~~~~~~~~~~~~~~#
#~~ functions ~~#
#~~~~~~~~~~~~~~~#

n_closest <- function(window,poly,pt,mt){
	n = nrow(mt)
	temp = matrix(seq(1:n),n,3) #1st col distance, 2nd col reference
	temp[,2] = mt[,1] #the x values
	temp[,3] = mt[,2] #the y values
	for(i in 1:n){
		temp[i,1] = abs(pt-mt[i,1]) #absolute value of points
	}
	temp = temp[order(temp[,1]),]
	temp = temp[1:window,]
	temp[,1] = temp[,1]/temp[window,1] #the max
	temp[,1] = tricube_w(temp[,1])
	W = diag(temp[,1])
	X = matrix(1,window,1)
	for(j in 1:poly){
		X = cbind(X,temp[,2]^poly)
	}
	Y = matrix(temp[,3],window,1)
	B = solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%Y
	Predicted = B[1,1] + B[2,1]*pt
	return(Predicted)
}

tricube_w <- function(tempVector){
	n = length(tempVector)
	temp = seq(1:n)
	for(i in 1:n){
		if (abs(tempVector[i]<1)){
			temp[i] = (1-(abs(tempVector[i]))^3)^3
		} else {
			temp[i] = 0
		}
	}
	return(temp)
}

distancex <- function(matrix,ref){
	n = nrow(matrix)
	temp = seq(1:n)
	for(i in 1:n){
		temp[i] = abs(-1*matrix[ref,1] + matrix[i,1])
	}
	return(temp)
}

F_scaledDistance <- function(x,s,n){
	temp_max = max(x[s:(s+n-1)])
	temp_array = x[s:(s+n-1)]/temp_max
	return(temp_array)
}

F_rows <- function(T_array,position,window_length){
	#this function finds the set of Xs that belong in the local set
	#local set has window_length elements in it
	left_n = position - window_length + 1
	right_n = position + window_length - 1
	n = length(T_array)
	if(left_n <= 0){ left_n = 1}
	if(right_n > n){right_n = n}
	new_m = matrix(0,right_n - left_n + 1,2)
	new_m[,1] = seq(left_n,right_n,by=1)
	new_m[,2] = T_array[left_n:right_n]
	new_m = new_m[order(new_m[,2]),]
	new_m2 = new_m[1:window_length,]
	new_m2 = new_m2[order(new_m2[,1]),]
	return(new_m2)
}

F_WLS <- function(matrixT,regressors, w, s, window)
{
	#weighted least squares. If F_lowess crashes, it's probably here (singular matrix)
	#this is linear. some fit higher order polynomial functions
	Y = matrix(matrixT[,2],window,1)
	X = regressors
	W = diag(w)
	B = solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%Y
	#Predicted = X%*%B
	return(B)
}

F_lowess <- function(d_data,window,iter=1,poly=1,density=FALSE){
	#d_data = matrix with x on 1st col, y on 2nd col
	#window = local area size
	#density=FALSE by default, make lowess denser  with more X values
	#d_data should be sorted by x
	lowess = d_data
	if(density==FALSE){
	for(jj in 1:iter){
		if(jj > 1){d_data = lowess}		
		n = nrow(d_data)
		for(ii in 1:n){
			position = ii
			distanceFromX = distancex(d_data,position)
			temp1 = F_rows(distanceFromX,position,window)
			localSet = d_data[temp1[,1],]
			X = matrix(0,window,1)
			X[,1] = 1
			for(kk in 1:poly){
				X = cbind(X,localSet[,1]^kk)
			}
			scaledDistance = F_scaledDistance(distanceFromX,temp1[1,1],window)
			weights = tricube_w(scaledDistance)
			wls_betas = F_WLS(localSet,X,weights,position,window) #if something goes wrong it's here
			#predicted = as.matrix(localSet)%*%wls_betas
			predicted = X%*%wls_betas
			rowN = which(temp1[,1]==position)
			lowess[ii,2] = predicted[rowN]
		}
	}
		return(lowess)
	}else{
		
		#7/24/2020 update
		#this makes it denser, none of this is optimized for speed
		nlength = nrow(d_data)
		dense = 4
		min_pt = min(d_data[,1])
		max_pt = max(d_data[,1])		
		denser = seq(min_pt,max_pt,(max_pt-min_pt)/(nlength*dense))
		denser = denser[2:(length(denser)-1)]
		lowess = matrix(-99999,length(denser),2)		
		lowess[,1] = denser		
		lowess = rbind(as.matrix(d_data),lowess)
		lowess = lowess[order(lowess[,1]),] #orders it by x values
		nn = nrow(lowess)
		#print(lowess)
		for(ii in 1:nn){
			lowess[ii,2] = n_closest(window,poly,lowess[ii,1],d_data)
			#print(paste("IN HERE",ii))	
		}
		return(lowess)
	}
}

#~~~~~~~~~~~~~~~~~#
#	Example	#
#~~~~~~~~~~~~~~~~~#

d_data = read.csv("E:/D/data/lowess.csv") #change file here for NIST data csv
select = 1 #1 NIST data, 2 for Sunspots, 3 for cars
zz = matrix(sunspot.month,length(sunspot.month),2)
zz[,1] = seq((1749+(1/12)),2013.75,by=(1/12))

#NIST DATA
if(select==1){
window = 7
poly = 1
it =1
lowess1 = F_lowess(d_data,window,it,poly,FALSE)
lowess2 = F_lowess(d_data,window,it,poly,TRUE)
title1 = paste("Example of LOWESS Smoother from nist.gov, iteration = ",it,", polynomial order =",poly)
plot(d_data[,1],d_data[,2],pch=19,main=title1, xlab="x",ylab="y")
lines(lowess1[,1],lowess1[,2],col="orange",lty=1,lwd=2)
lines(lowess2[,1],lowess2[,2],col="firebrick",lty=1,lwd=2)
}

#Sunspots monthly
if(select==2){
d_data = zz
window = 12#in months
window2 = 120*4
it = 3 #iterations
it2= 1
poly = 1
lowess = F_lowess(d_data,window,it,poly,FALSE)
lowess2 = F_lowess(d_data,window2,it2,poly,FALSE)
lowess3 = F_lowess(d_data[1:2540,],window2,it2,poly,FALSE)
title1 = paste("Sunspots with LOWESS smoother" )
plot(d_data[,1],d_data[,2],type="p",pch=19,main=title1,xlab="X (years)",ylab="Y (sunspots)",col="black")
lines(d_data[,1],lowess[,2],col="orange",lty=1,lwd=3)
lines(d_data[,1],lowess2[,2],col="firebrick",lty=1,lwd=3)
lines(lowess3[,1],lowess3[,2],col="firebrick",lty=2,lwd=3)
legend1 = paste("window =",window,"months, iteration =",it)
legend2 = paste("window =",window2,"months, iteration =",it2)
#legend("topleft",bty="n",c(legend2,legend2),lty=c(1,2),col=c("firebrick","firebrick"),lwd=3)
legend("topleft",bty="n",c(legend1,legend2),lty=1,col=c("orange","firebrick"),lwd=3)
}


if(select==3){
zz3 = matrix(cars$speed,50,2)
zz3[,2] = cars$dist
poly = 1
cars1 = F_lowess(zz3,20,1,poly,TRUE)
cars2 = F_lowess(zz3,20,2,poly,FALSE)
plot(cars$speed, cars$dist, main="cars with window = 20",pch=19)
lines(cars1[,1],cars1[,2],col="blue",lwd=2)
lines(cars2[,1],cars2[,2],col="purple",lwd=2)
legend("topleft", c(paste("iteration =", c("1", "2"))), lty=1,lwd=2, col=c("blue","purple"),bty="n")
}


