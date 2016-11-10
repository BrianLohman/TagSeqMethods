cat("Plotting functions:\nstrip -- show the average as a point with option SE / CI, shows the data\nbar -- barplot, still needs some work\nscatter -- scatterplot\nmod -- modified box plot\nhistogram -- self-explanatory\nsimple -- draws a line for the mean of each group and plots the data points\nbeeStrip -- like simple, but plotted using the beeswarm package so no alpha is needed\naddAlpha -- add transparancy to any color. pass the color and alpha value in [0,1]\n\nStats functions:\nmonte_unpaired -- monte carlo permutation test, two unpaired samples\nmonte_paired -- monte carlo permutation test, two paired samples\n")  
if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
if (!"wesanderson" %in% installed.packages()) install.packages("wesanderson")
if (!"magrittr" %in% installed.packages()) install.packages("magrittr")
if (!"beeswarm" %in% installed.packages()) install.packages("beeswarm")
library(beeswarm)
library(dplyr)
library(wesanderson)
library(magrittr)
library(viridis)

# some attractive colors:
ruby <-rgb(202/255,53/255,7/255,1)
slate <- rgb(137/255,157/255,164/255,1)
mint <- rgb(73/255,191/255,150/255,1)
golden <- rgb(218/255,165/255,32/255,1)
orange <- rgb(227/255,121/255,46/255,1)
sky <- rgb(95/255,188/255,212/255,1)

cols<-c(ruby,mint,golden,slate,orange,sky)



# things to look into:
# alluvial diagrams: http://www.r-bloggers.com/alluvial-diagrams/
# parcoord in MASS
# stars()


## line for mean + data
## inpout data as a list; a mean is returned for each element of list. add median == T to add a diamond median
## example:
# x <- list(rnorm(40,40,5),rnorm(20,35,2),rnorm(25,41,2))
# simple(x,main="simple() defaults") # using the defaults
### simple(x,jitter=F) # without jitter doesn't look as good
### simple(x,line_col="black",point_col=c(ruby,mint,slate),ylab="measurement",xlab="group",lab=c("A","B","C"),rug=T)
## add the IQR: simple(x,line_col="black",point_col=c(ruby,mint,slate),ylab="measurement",xlab="group",lab=c("A","B","C"),rug=T,IQR=T)
# simple(list(rnorm(50,50,5),rnorm(30,40,6),rnorm(10,60,2),rnorm(60,50,10),rnorm(30,39,4)),point_col=wes_palette(5, name = "Zissou"),line_color="black",median=T,main="simple() dressed up")
## simple(list(rnorm(24,10,2),rchisq(20,5),rexp(40,1/5),runif(40,5,15)),lab=c("normal","chi-squared","exponetial","uniform"),point_col=c(ruby,mint,slate,"goldenrod"),line_color="black",median=T)

# to do: replace rug place with marginal density plots of each group?

simple<-function(data,lab=rep(c(),length(data)),point_size=1.2,line_color="red",line_width=3.0,jitter=T,point_col=viridis(length(data)+2)[1:length(data)],y_limits=c(min(unlist(data),na.rm=T),max(unlist(data),na.rm=T)),median=FALSE,rug=TRUE,sample_size=T,IQR=F,...){
  
  # make sure the data is a list; in the future I need to change the code to accomodate other types of data input (e.g. dataframes with y~x)
  if(is.list(data)!=TRUE){
    stop("input your data as a list")
  }
  
  number_groups<-length(data)
  
  # to get nice x values for plotting, we'll use those generate by barplot()
  x_values<-barplot(rep(1,number_groups),plot=F)
  
  # create the plot
  plot(c(x_values[1]*0.4,x_values[length(x_values)]*1.2),y_limits,type="n",xaxt="n",bty="l",...)
  
  # find range of x:
  if(length(data)==1){
    xRange <- 0.5
  }
  else{
    xRange <- range(x_values,na.rm=T)[2]-range(x_values,na.rm=T)[1]
  }
  
  # plot the mean as a line
  for(i in 1:length(data)){
    lines(x=c(x_values[i]-(0.1*xRange),x_values[i]+(0.1*xRange)),y=c(rep(mean(data[[i]],na.rm=T),2)),col=line_color,lwd=line_width)
  }
  
  # to plot the IQR:
  if(IQR==TRUE){
    p<-boxplot(data,plot=F)
    for(i in 1:length(data)){
      lines(x=c(x_values[i],x_values[i]),y=c(p$stats[2,i],p$stats[4,i]),lwd=2,col=addAlpha(line_color))
    }
  }
  
  # create x axis
  axis(side=1,at=x_values,labels=lab)
  
  # make point_col length equal to number of groups
  if(length(point_col==1)){
    point_col = rep(point_col,number_groups)
  }
  
  # draw data points 
  if(jitter==T){
    for(i in 1:number_groups){
      points(x=rep(x_values[i],length(data[[i]]))+rnorm(length(data[[i]]),0,xRange*0.02),y=data[[i]],pch=16,col=point_col[i],cex=point_size)
    }
  }
  else{
    for(i in 1:number_groups){
      points(x=rep(x_values[i],length(data[[i]])),y=data[[i]],pch=16,col=point_col[i],cex=point_size)
    }
  }
  
  # for rug plotting:
  if(rug==T){
    for(i in 1:length(data)){
      rug(data[[i]],side=4,col=addAlpha(point_col[i],0.4),lwd=2)
    }
  }
  
  # for plotting the median as a diamond:
  if(median==T){
    points(x=x_values,y=unlist(lapply(data,median)),pch=5,col=line_color,cex=1.5,lwd=3)
  }
  
  # for plotting sample size below each group
  if(sample_size==TRUE){
    for(i in 1:length(data)){
      text(x_values[i],par("usr")[3]*1.02,paste("n = ",length(data[[i]]),sep=""),col="grey20")
    }
  }
}

# similar to simple() above, but plot using beeswarm

# examples
## beeStrip(list(rnorm(50,50,5),rnorm(30,40,6),rnorm(10,60,2),rnorm(60,50,10),rnorm(30,39,4)),point_col=wes_palette(5, name = "Zissou"),line_color="black",median=T)
# beeStrip(list(iris %>% filter(Species=="setosa") %>% .$Sepal.Length, iris %>% filter(Species=="versicolor") %>% .$Sepal.Length, iris %>% filter(Species=="virginica") %>% .$Sepal.Length),point_col=viridis(3),line_color="black",IQR=T,lab=c("setosa","versicolor","virginica"),xlab="species",ylab="sepal length")
# beeStrip(list(rnorm(50,50,5),rnorm(30,40,6),rnorm(10,60,2),rnorm(60,50,10),rnorm(30,39,4)),point_col=wes_palette(5, name = "Zissou"),line_color="grey30",IQR=T)

beeStrip<-function(data,lab=rep(c(),length(data)),point_size=1.4,beeMethod="center",line_color="red",line_width=3.0,jitter=T,point_col=viridis(length(data)+1)[1:length(data)],y_limits=c(min(unlist(data),na.rm=T),max(unlist(data),na.rm=T)),median=FALSE,rug=TRUE,sample_size=T,IQR=F,side=-1,...){
  
  # make sure the data is a list; in the future I need to change the code to accomodate other types of data input (e.g. dataframes with y~x)
  if(is.list(data)!=TRUE){
    stop("input your data as a list")
  }
  
  number_groups<-length(data)
  
  # create the plot
  beeswarm(data,method=beeMethod,priority="density",pch=16,col=point_col,cex=point_size,side = side,bty='l',labels=lab,...)
  
  x_values = 1:number_groups
  
  # find range of x:
  xRange <- length(data)
  
  # plot the mean as a line
  for(i in 1:length(data)){
    scale = 0.05
    lines(x=c(x_values[i]-(scale*xRange),x_values[i]+(scale*xRange)),y=c(rep(mean(data[[i]],na.rm=T),2)),col=line_color,lwd=line_width)
  }
  
  # to plot the IQR:
  if(IQR==TRUE){
    p<-boxplot(data,plot=F)
    for(i in 1:length(data)){
      lines(x=c(x_values[i],x_values[i]),y=c(p$stats[2,i],p$stats[4,i]),col=line_color,lwd=line_width)
    }
  }
  
  # create x axis
  axis(side=1,at=x_values,labels=lab)
  
  # make point_col length equal to number of groups
  if(length(point_col==1)){
    point_col = rep(point_col,number_groups)
  }
  
  # for rug plotting:
  if(rug==T){
    for(i in 1:length(data)){
      rug(data[[i]],side=4,col=addAlpha(point_col[i],0.4),lwd=4)
    }
  }
  
  # for plotting the median as a diamond:
  if(median==T){
    points(x=x_values,y=unlist(lapply(data,median)),pch=5,col=line_color,cex=1.5,lwd=3)
  }
  
  # for plotting sample size below each group
  if(sample_size==TRUE){
    for(i in 1:length(data)){
      text(x_values[i],y_limits[1]*0.99,paste("n = ",length(data[[i]]),sep=""),col="grey20")
    }
  }
}



#data should be a list of vectors you want to plot

# examples
# strip(list(rnorm(50,5),rnorm(20,10)))
# strip(list(rpois(50,5),rpois(20,3),rnorm(100,7)),type="ci",lab=c("a","b","c"),xlab="group",ylab="response")
# lots of groups:
#strip(list(rpois(50,5),rpois(20,3),rnorm(100,7),rnorm(45,2),rnorm(90,9)),type="ci",lab=c("a","b","c","d","e"),xlab="group",ylab="response")

strip<-function(data,lab=rep(c(),length(data)),type="se",jitter=T,points=16,xlab="",ymin="determine",ymax="determine",point_size=1.2,mean_col = "red",cols = viridis(length(data)+2)[1:length(data)],...){
  par(bty="l")
  
  # error checking:
  # if there are NAs in your dataset, throw an error
  if(length(Filter(is.na,data))>0){
    stop("there are NAs in your list")
  }
  
  number_groups<-length(data)
  factor_names<-lab
  
  
  if(ymin=="determine"&ymax=="determine"){
    maximum_value = -1000
    minimum_value = 1000
    for(i in 1:length(data)){
      if(max(data[[i]]>maximum_value)){
        maximum_value <- max(data[[i]])
      }
      if(min(data[[i]]<minimum_value)){
        minimum_value <- min(data[[i]])
      }
    }
  }
  
  else{
    maximum_value=ymax
    minimum_value=ymin
  }

  
  plot(c(0,1),c(minimum_value,maximum_value),type="n",xaxt="n",xlab=xlab,...)
  
  if(number_groups == 1){
    x_values = 0.5
    offset= 0.15
  }
  else{
    x_values<-seq(0.2,0.8,length.out=number_groups)
    # offset is the distance between the data points and their means
    offset <- 0.2 / number_groups
  }
  
  # use lapply to extract the means
  means<-lapply(data,mean)
  
  
  if(type=="se"){
    std_dev<-lapply(data,sd)
    for(i in 1:length(data)){
      std_error <- std_dev[[i]] / sqrt(length(data[[i]]))
      arrows(x_values[i]+offset,means[[i]]-std_error,x_values[i]+offset,means[[i]]+std_error,angle=90,code=3,length=0.07)
    }
  }

  if(type=="sd"){
    std_dev<-lapply(data,sd)
    for(i in 1:length(data)){
      arrows(x_values[i]+offset,means[[i]]-std_dev[[i]],x_values[i]+offset,means[[i]]+std_dev[[i]],angle=90,code=3,length=0.07)
    }
  }
  
  if(type=="ci"){
    std_dev<-lapply(data,sd)
    for(i in 1:length(data)){
      std_error <- std_dev[[i]] / sqrt(length(data[[i]]))
      arrows(x_values[i]+offset,t.test(data[[i]])$conf[1],x_values[i]+offset,t.test(data[[i]])$conf[2],angle=90,code=3,length=0.07)
    }
  }
  
  # plot the means
  for(i in 1:length(data)){
    points(x=x_values[i]+offset,y=means[[i]],cex=point_size*1.3,pch=16,col=mean_col)
  }
  
  
  axis(side=1,at=x_values,labels=lab)
  
  
  if(jitter==F){
    # now to draw the points
    for(i in 1:number_groups){
      points(x=rep(x_values[i],length(data[[i]])),y=data[[i]],pch=points,col=cols[i])
    }
  }
  
  else if(jitter==T){
    for(i in 1:number_groups){
      points(x=rep(x_values[i],length(data[[i]]))+rnorm(length(data[[i]]),0,0.02/(number_groups/2)),y=data[[i]],pch=points,col=cols[i],cex=point_size)
    }
  }
  par(bty="o",lwd=1)
}


#########################################################


# strip2 works fine, but it only accepts dataframes
# to increase flexibility, I've introduced strip(), which takes a list of vectors you want to plot
# it's a lot more flexible
strip2<-function(response,group,data=dat,lab=levels(group),type="se",jitter=T,points=16,...){
  par(bty="l",lwd=1.7)
  if(is.factor(group)!=TRUE){
    warning("'group' is not a factor: converting now")
    group<-factor(group)
  }
  
  number_groups<-nlevels(group)
  factor_names<-levels(group)
  boxplot_table<-boxplot(response~group,plot=F)
  
  # colors:
  cols<-c("#3B9AB2","#F21A00","#5B1A18","#E1AF00","#446455")
  cols_points<-c("#3B9AB270","#F21A0070","#5B1A1870","#E1AF0070","#44645570")
  
  plot(c(0,1),c(min(response),max(response)),type="n",xaxt="n",...)
  
  if(number_groups == 1){
    x_values = 0.5
    offset= 0.15
  }
  else{
    x_values<-seq(0.2,0.8,length.out=number_groups)
    # offset is the distance between the data points and their means
    offset <- 0.15 / number_groups
  }
  

  
  # use tapply to extract the means
  means<-tapply(response,group,mean)
  
  
  if(type=="se"){
    # now for +- one s.e.
    std_dev<-tapply(response,group,sd)
    std_error <- std_dev / sqrt(boxplot_table$n)
    # draw standard errors:
    for(i in 1:number_groups){
      arrows(x_values[i]+offset,means[i]-std_error[i],x_values[i]+offset,means[i]+std_error[i],angle=90,code=3,length=0.07)
    }
  }
  
  if(type=="sd"){
    # now for +- one s.e.
    std_dev<-tapply(response,group,sd)
    # draw standard deviations:
    for(i in 1:number_groups){
      arrows(x_values[i]+offset,means[i]-std_dev[i],x_values[i]+offset,means[i]+std_dev[i],angle=90,code=3,length=0.07)
    }
  }
  
  if(type=="ci"){
    std_dev<-tapply(response,group,sd)
    std_error <- std_dev / sqrt(boxplot_table$n)
    for(i in 1:number_groups){
      arrows(x_values[i]+offset,means[i]-(1.96*std_error[i]),x_values[i]+offset,means[i]+(1.96*std_error[i]),angle=90,code=3,length=0.07)
    }
  }
  
  if(type=="cred"){
    cred<-tapply(response,group,credible_intervals)
    print("error bars are credible intervals determined by 100,000 permutations")
    for(i in 1:number_groups){
      arrows(x_values[i]+offset,cred[[i]][[1]],x_values[i]+offset,cred[[i]][[2]],angle=90,code=3,length=0.07)
    }
  }
  
  # plot the means
  points(x=x_values+offset, y=means,cex=1.5,pch=16,col=cols)
  
  
  axis(side=1,at=x_values,labels=lab)
  
  
  if(jitter==F){
    # now to draw the points
    for(i in 1:number_groups){
      points(x=rep(x_values[i],boxplot_table$n[i]),y=response[group==factor_names[i]],pch=points,col=cols_points[i])
    }
  }
  
  else if(jitter==T){
    for(i in 1:number_groups){
      points(x=rep(x_values[i],boxplot_table$n[i])+rnorm(boxplot_table$n[i],0,0.02/(number_groups/2)),y=response[group==factor_names[i]],pch=points,col=cols_points[i])
    }
  }
  par(bty="o",lwd=1)
}


# for making credible intervals
# used in strip() above
credible_intervals<-function(x,n=100000){
  means<-vector(length=n)
  for(i in 1:n){
    means[i]<-mean(sample(x,size=length(x),replace=T))
  }
  #print("first result is 5th percentile, second is 95th")
  results<-quantile(means,c(0.05,0.95))
  print(results)
  #hist(means)
  return(results)
}

##############################################

# barplot

# sim is a list of vectors you want to plot
# lab is the labels for under each bar

# examples of usage:
# simplest:
# bar(list(rnorm(20,5),rnorm(15,5)),lab=c("A","B"),xlab="group",ylab="response",jitter=F)
# with conf. intervals:
# bar(list(rnorm(20,5),rnorm(15,5)),CI=T,lab=c("A","B"),xlab="group",ylab="response",jitter=F)
# nicer:
# bar(list(rnorm(20,5),rnorm(15,5)),CI=T,lab=c("A","B"),xlab="group",ylab="response")
# with more groups
# bar(list(rnorm(20,5),rnorm(15,5),rnorm(200,50),rnorm(50,1),rnorm(34,19)),CI=T)
# with medians plotted:
# bar(list(rnorm(20,5),rnorm(15,5),rnorm(200,50),rnorm(50,1),rnorm(34,19)),CI=T,median=T)
bar<-function(sim,lab=rep(c(),length(sim)),CI=F,SE=F,bar_color="grey80",jitter=T,point_col="#00000080",y_limits=NA,median=FALSE,...){
  par(lwd = 1,family = 'Helvetica')
  
  if(is.list(sim)==F){
    stop("the first argument needs to be a list of the vectors you want to plot")
  }
  
  means<-vector(length=length(sim))
  for(i in 1:length(sim)){
    means[i]<-mean(sim[[i]])
  }
  
  max_value = 0
  for(i in 1:length(sim)){
    if(max(sim[[i]],na.rm=TRUE) >= max_value){
      max_value <- max(sim[[i]]) 
    }
  }
  
  min_value = 0
  for(i in 1:length(sim)){
    if(min(sim[[i]],na.rm=TRUE) <= min_value){
      min_value <- min(sim[[i]]) 
    }
  }
  
  if(is.na(y_limits)==T){
    (p<-barplot(means,bty="l",ylim=c(floor(min_value),ceiling(max_value)+(0.1*max_value)),axes=F,col=bar_color,border=NA,...))
  }
  else{
    (p<-barplot(means,bty="l",ylim=y_limits,axes=F,col=bar_color,border=NA,...))
  }
  #(p<-barplot(means,bty="l",space=0.4,ylim=ifelse(y_limits==NA,c(floor(min_value),ceiling(max_value)+(0.1*max_value)),y_limits),axes=F,col=bar_color,border=NA,...))
  axis(1,at=p,lwd=1,cex=1.5,labels=lab,tick=F)
  axis(2,cex=1.5,lwd=1)
  
  if(SE==T & CI ==T){
    stop("can't plot both a conf. interval and a standard error. choose only one.")
  }
  
  # for 95% CIs:
  # why doesn't this work?
  if(CI==T){
    for(i in 1:length(sim)){
      # find the stardard error
      arrows(x0=p[i],y0=t.test(sim[[i]])$conf[1],x1=p[i],y1=t.test(sim[[i]])$conf[2],col="grey50",angle=90,code=3,lwd=1,length=.1)
      #lines(rep(p[i],2),y=c(mean(sim[[i]]),mean(sim[[i]])-(se*1.96)),lwd=5,col="white")
      #lines(rep(p[i],2),y=c(mean(sim[[i]]),mean(sim[[i]])+(se*1.96)),lwd=5)
    }
  }
  
  if(SE==T){
    for(i in 1:length(sim)){
      # find the stardard error
      se<-sd(sim[[i]])/sqrt(length(sim[[i]]))
      arrows(x0=p[i],y0=mean(sim[[i]])-se,x1=p[i],y1=mean(sim[[i]])+se,col="grey50",angle=90,code=3,lwd=1,length=.1)
      #lines(rep(p[i],2),y=c(mean(sim[[i]]),mean(sim[[i]])-se),lwd=5,col="white")
      #lines(rep(p[i],2),y=c(mean(sim[[i]]),mean(sim[[i]])+se),lwd=5)
    }
  }
  
  if(jitter==T){
    x_vals = p + 0.2
    for(i in 1:length(sim)){
      points(x=rep(x_vals[i],length(sim[[i]]))+rnorm(length(sim[[i]]),0,0.05),y=sim[[i]],pch=16,col=point_col)
    }
  }
  
  # plot median plots
  if(median==T){
    medians<-vector(length=length(sim))
    for(i in 1:length(sim)){
      medians[i]<-median(sim[[i]],na.rm=T)
    }
    points(p,medians,pch=5)
  }
  
  return(p)
}

########################################################
#scatter


# helper function that will draw line of best fit, but not extend beyond the reaches of the data
line<-function(x,y,color="black",...){
  slope <- summary(lm(y~x))$coefficients[2,1]
  intercept <- summary(lm(y~x))$coefficients[1,1]
  lines(x=c(min(x),max(x)),y=c(transform(min(x),m=slope,b=intercept),transform(max(x),m=slope,b=intercept)),col=color,...)
}

# returns y given x along a line
transform <- function(x,m,b){
  return((m*x)+b)
}


scatterOld<-function(x,y,xlab="",ylab="",line=T,sig=T,color="black",line_col="red",confidenceInterval=T,...){
  par(lwd=2,cex=1,bg="white")
  # error checking
  if(length(x)!=length(y)){
    stop("x and y lengths differ")
  }
  
  # do the actual plotting
  plot(x,y,xlab=xlab,ylab=ylab,pch=16,bty="n",col=color,...)
  
  p<-summary(lm(y~x))$coefficients[2,4]
  c<-summary(lm(y~x))$coefficients[2,1]
  
  # draw the line of best fit
  if(line==T){
    if(p<=0.05){
      line(x,y,lwd=2,color=line_col)
    }
    else{
      line(x,y,lwd=2,lty=2,color=line_col)
    }
  }
  
  # add statistical significance to the plot
  if(sig==T){
    if(c<0){
      print(c)
      legend(x="bottomright",bty="n",legend=paste("p = ",round(p,3)))
    }
    else{
      print(c)
      legend(x="bottomleft",bty="n",legend=paste("p = ",round(p,3)))
    }
  }
  
  # adding confidence intervals around the regression line
  if(confidenceInterval==T){
    model<-lm(y~x)
    xVals<-seq(min(x),max(x),.1)
    conf<-predict(model,data.frame(x=xVals),interval="confidence",level=0.95)
    lines(xVals,conf[,2],col="#00000050",lty=2)
    lines(xVals,conf[,3],col="#00000050",lty=2)
  }
  
  # reset the par values
  par(lwd=1,cex=1,bg="white")
}



# example
## scatter(rnorm(50,50,10),rnorm(50,30,3))

scatter<-function(x,y,xlab="",ylab="",line=T,sig=T,color="black",line_col="red",confidenceInterval=T,plottingCharacter=16,rug = T,...){
  par(lwd=1,cex=1,bg="white")
  
  # error checking
  if(length(x)!=length(y)){
    stop("x and y lengths differ")
  }
  
  # do the actual plotting
  plot(x,y,xlab=xlab,ylab=ylab,pch=plottingCharacter,bty="l",col=color,...)
  
  p<-summary(lm(y~x))$coefficients[2,4]
  c<-summary(lm(y~x))$coefficients[2,1]
  
  # draw the line of best fit
  if(line==T){
    if(p<=0.05){
      line(x,y,lwd=2,color=line_col)
    }
    else{
      line(x,y,lwd=2,lty=2,color=line_col)
    }
  }
  
  # adding confidence intervals around the regression line
  if(confidenceInterval==T){
    model<-lm(y~x)
    xVals<-seq(min(x),max(x),.1)
    conf<-predict(model,data.frame(x=xVals),interval="confidence",level=0.95)
    lines(xVals,conf[,2],col="#00000050",lty=2,lwd=2)
    lines(xVals,conf[,3],col="#00000050",lty=2,lwd=2)
  }
  
  # adding 1-d 'rug plots' to the bottom and right of the plot to show distribution of each variable
  if(rug==T){
    rug(x,side=1,col="#00000070",lwd=2)
    rug(y,side=2,col="#00000070",lwd=2)
  }
  
  # add statistical significance to the plot
  if(sig==T){
    if(c>0){
      legend(x="bottomright",bty="n",legend=paste("p = ",round(p,3)))
    }
    else{
      legend(x="bottomleft",bty="n",legend=paste("p = ",round(p,3)))
    }
  }
  
  # reset the par values
  par(lwd=1,cex=1,bg="white")
}

##############
### beeStrip + mod
##################
#data(iris)
#beeStripMod(iris$Sepal.Length,iris$Species,xlab="species",ylab="sepal length",main="beeStripMod() example")

beeStripMod<-function(data,group,lab=rep(c(),length(data)),point_size=1.4,beeMethod="center",line_width=3.0,jitter=T,point_col=ifelse(is.list(data) %>% rep(20),viridis(length(data)+1)[1:length(data)],viridis(nlevels(group)+1)[1:nlevels(group)]),y_limits=c(ifelse(is.list(data),min(unlist(data)),min(data)),ifelse(is.list(data),max(unlist(data),max(data)))),mean=FALSE,sample_size=T,side=-1,red_median=F,stats=T,...){
  
  # if response is missing, assume data is a list
  if(missing(group)){
    if(is.list(data)==FALSE){
      stop("enter your data either as a list or as a response variable and factor within a dataframe")
    }
  }
  
  if(is.list(data)){
    print("data are a list")
    number_groups<-length(group)
    # get statistics you need to plot
    boxplot_table<-boxplot(data,plot=F)
    # create the plot
    beeswarm(data,method=beeMethod,priority="density",pch=16,col=point_col,cex=point_size,side = side,bty='l',labels=lab,...)
  }
  
  else{
    number_groups<-nlevels(group)
    boxplot_table<-boxplot(data~group,plot=F)
    # create the plot
    beeswarm(data~group,method=beeMethod,priority="density",pch=16,col=point_col,cex=point_size,side = side,bty='l',labels=lab,...)
    
    if(stats==T){
      print(TukeyHSD(aov(data~group)))
      x = summary(aov(data~group))[[1]][["Pr(>F)"]][1] %>% round(4)
      if(x ==0){
        x = "p < 0.001"
      }
      else{
        x = paste("p = ",x,sep="")
      }
      text(x=(median(1:number_groups)),y=max(data),labels=paste("anova, ",x,sep=""),pos=1)
    }
    
  }
  
  x_values = 1:number_groups
  point_col = point_col[1:number_groups]
  scale = 0.05
  
  # plot the mean as a point
  for(i in 1:number_groups){
    points(x=x_values[i]+(scale*number_groups), y=boxplot_table$stats[3,i],cex=1.1,pch=16,col=ifelse(red_median==T,"red",point_col[i]))
  }
  # draw lines representing 1 and fourth quartiles
  for(i in 1:number_groups){
    lines(x=rep(x_values[i]+(scale*number_groups),2),y=c(boxplot_table$stats[1,i],boxplot_table$stats[2,i]),lwd=2)
    lines(x=rep(x_values[i]+(scale*number_groups),2),y=c(boxplot_table$stats[4,i],boxplot_table$stats[5,i]),lwd=2)
  }
  
  
  # for plotting sample size below each group
  if(sample_size==TRUE){
      text(x=x_values+(scale*number_groups),y=min(data)*0.99,paste("n = ",boxplot_table$n,sep=""),col="grey20")
  }
}


##############
### beeStrip + boxplot
##################
#data(iris)
#beeStripBox(iris$Sepal.Length,iris$Species,xlab="species",ylab="sepal length",main="beeStripBox() example")
# TO DO: fix box color argument. can't get boxes to take on the right colors
beeStripBox<-function(data,group,lab=rep(c(),length(data)),point_size=1.4,beeMethod="center",line_width=3.0,point_col=ifelse(is.list(data) %>% rep(20),viridis(length(data)+1)[1:length(data)],viridis(nlevels(group)+1)[1:nlevels(group)]),y_limits=c(ifelse(is.list(data),min(unlist(data)),min(data)),ifelse(is.list(data),max(unlist(data),max(data)))),mean=FALSE,sample_size=T,side=-1,red_median=F,stats=T,box_thickness = 0.2,box_color=FALSE,...){
  
  # if response is missing, assume data is a list
  if(missing(group)){
    if(is.list(data)==FALSE){
      stop("enter your data either as a list or as a response variable and factor within a dataframe")
    }
  }
  
  if(is.list(data)){
    print("data are a list")
    number_groups<-length(data)
    point_col = point_col[1:number_groups]
    # get statistics you need to plot
    boxplot_table<-boxplot(data,plot=F)
    # create the plot
    beeswarm(data,method=beeMethod,priority="density",pch=16,col=point_col,cex=point_size,side = side,bty='l',labels=lab,...)
    boxplot(data, main = "", axes = FALSE, at = 1:number_groups+0.2, xlab=" ", ylab=" ", border = ifelse(box_color==T,point_col,"black"), add=TRUE ,boxwex = box_thickness,pars = list(medlty = 1, whisklty = c(1, 1), medcex = 0.7, outcex = 0, staplelty = "blank"))
  }
  
  else{
    number_groups<-nlevels(group)
    boxplot_table<-boxplot(data~group,plot=F)
    point_col = point_col[1:number_groups]
    # create the plot
    beeswarm(data~group,method=beeMethod,priority="density",pch=16,col=point_col,cex=point_size,side = side,bty='l',labels=lab,...)
    boxplot(data~group, add=TRUE, main = "", at = 1:number_groups+0.2, axes = FALSE, xlab=" ", ylab=" ",border = ifelse(box_color==T,point_col,"black"), boxwex = box_thickness, pars = list(medlty = 1, whisklty = c(1, 1), medcex = 0.7, outcex = 0, staplelty = "blank"))
    
    if(stats==T){
      print(TukeyHSD(aov(data~group)))
      x = summary(aov(data~group))[[1]][["Pr(>F)"]][1] %>% round(4)
      if(x ==0){
        x = "p < 0.001"
      }
      else{
        x = paste("p = ",x,sep="")
      }
      text(x=(median(1:number_groups)),y=max(data),labels=paste("anova, ",x,sep=""),pos=1)
    }
    
  }
  
  x_values = 1:number_groups
  scale = 0.05
  
  # for plotting sample size below each group
  if(sample_size==TRUE){
    text(x=x_values+(scale*number_groups),y=min(data)*0.99,paste("n = ",boxplot_table$n,sep=""),col="grey20")
  }
}




## modified boxplot
## example:
### x<-data.frame(c(rnorm(20,5),rnorm(30,10,1.5),rnorm(25,12,1.6),rnorm(20,20,3)),c(rep("A",20),rep("B",30),rep("C",25),rep("D",20)))
### mod(x[,1],x[,2])
mod<-function(response,group,lab=levels(group),...){
  
  # make sure group is a factor; if not, convert it to one
  if(is.factor(group)!=TRUE){
    warning("'group' is not a factor: converting now")
    group<-factor(group)
  }
  
  number_groups<-nlevels(group)
  
  # get statistics you need to plot
  boxplot_table<-boxplot(response~group,plot=F)
  
  # create the plot
  plot(c(0,1),c(min(response),max(response)),type="n",xaxt="n",...)
  x_values<-seq(0.7,0.3,length.out=number_groups)
  
  # plot the median
  points(x=x_values, y=boxplot_table$stats[3,1:number_groups],cex=1.1,pch=16,col="red")
  
  # creare x axis
  axis(side=1,at=x_values,labels=lab)
  
  # draw lines representing 1 and fourth quartiles
  for(i in 1:number_groups){
    lines(x=rep(x_values[i],2),y=c(boxplot_table$stats[1,i],boxplot_table$stats[2,i]),lwd=2)
    lines(x=rep(x_values[i],2),y=c(boxplot_table$stats[4,i],boxplot_table$stats[5,i]),lwd=2)
  }
}

# histogram
histogram=function(x,color="grey50",bor="grey50",...){
  hist(x,col=color,border=bor,lwd=2,cex.lab=1.2,...)
}


### some non-parametric stats functions
#############################################
### monte carlo, two-samples, unpaired #####
###########################################
# the idea here is to randomize associations of the data with the group
# then calculate the difference in the group means under the null hypothesis
# then compare it to the observed difference in means

## should do some more error checking to make sure this works, but looks good, 14 may 2014


monte_unpaired<-function(group_a,group_b,n=9999,null=0,table=TRUE){
  values<-vector(length=n)
  crit = mean(group_a)-mean(group_b)
  diff<-abs(crit-null)
  
  # make a data frame
  dat<-data.frame(c(group_a,group_b),c(rep("group_a",length(group_a)),rep("group_b",length(group_b))))
  names(dat)<-c("values","group")
  
  for(i in 1:n){
    dat$group<-sample(dat$group,replace=F)
    x<-tapply(dat$value,dat$group,mean)
    values[i]<-x[[1]]-x[[2]]
  }
  
  lower_bound <- quantile(values,0.025)[[1]]
  upper_bound <- quantile(values,0.975)[[1]]
  
  p <- ((length(values[values<=null-diff])+length(values[values>=null+diff]))+1)/(n+1)
  # why add the +1's? see Ruxton and Neuha¨user, methods in ecology and evolution, 2013, doi:10.1111/2041-210X.12102
  
  # plot
  par(lwd = 3,family = 'Helvetica',cex.lab=1.3,cex.lab=1.3)
  (plot<-hist(values,cex.lab=1.3,xlab="Simulated Differences",main="Monte Carlo Simulation",col="#99999940",breaks=20))
  text(mean(values),max(plot$counts),paste("observed = ",round(crit,3)),pos=4)
  segments(crit,0,crit,n,lty=2)
  segments(lower_bound,0,lower_bound,n,lty=3,lwd=2)
  segments(upper_bound,0,upper_bound,n,lty=3,lwd=2)
  
  if(table==TRUE){
    g<-c(round(lower_bound,2),round(upper_bound,2),round(p,3))
    h<-c("lower_bound","upper_bound","p-value")
    result<-data.frame(h,g)
    colnames(result)<-c(" "," ")
    print(paste("based on ",n+1," iterations"))
    # I construct a confidence interval for the p-value based on the same paper above:
    # Ruxton and Neuha¨user, methods in ecology and evolution, 2013, doi:10.1111/2041-210X.12102
    ci <- 1.96*sqrt((p*(1-p))/(n+1))
    print(paste("lower confidence level for p-value:",round(p-ci,3),". Upper confidence interval: ",round(p+ci,3)))
    print(result)
    #return(values)
  }
  else{
    print(paste("based on",n,"iterations"))
    return(p) 
  }
}


#######################################################
# monte carlo method for a two-sample paired t-test ##
#####################################################
# example of use:
#monte_paired(df,n=9999)

monte_paired<-function(data,n=9999,table=TRUE,diff_col=3,one_sided=F){
  x<-vector(length=n)
  if(!is.data.frame(data)){
    stop("please enter a data frame as the first argument")
  }
  
  if(one_sided==T){
    crit <- (mean(data[,diff_col]))
  }
  else{
    crit <- abs(mean(data[,diff_col]))
  }
  
  
  for(i in 1:n){
    data1<-data.frame(t(apply(df[,1:2],1,sample)))
    data1$diff<-data1[,1]-data1[,2]
    if(one_sided==T){
      x[i]<-(mean(data1$diff))
    }
    else{
      x[i]<-abs(mean(data1$diff))
    }
  }
  
  
  extreme <- length(x[x>=crit])+1
  p <- extreme/(n+1)
  # why add the +1's? see Ruxton and Neuha¨user, methods in ecology and evolution, 2013, doi:10.1111/2041-210X.12102
  par(lwd = 3,family = 'Helvetica',cex.lab=1.3,cex.lab=1.3)
  
  (plot<-hist(x,cex.lab=1.3,xlab="Simulated Differences",main="MC Simulation",col="#99999940"))
  text(mean(x),max(plot$counts),paste("observed = ",round(crit,1)),pos=4)
  segments(crit,0,crit,n,lty=2)
  
  if(table==TRUE){
    g<-c(crit,round(extreme,0),p)
    h<-c("observed diff","num_extreme","p-value")
    result<-data.frame(h,g)
    colnames(result)<-c(" "," ")
    print(result)
    #return(x)
  }
  else{
    return(x) 
  }
}

## Add an alpha value to a colour
# from http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html
addAlpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}