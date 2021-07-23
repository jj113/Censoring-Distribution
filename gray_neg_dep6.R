
rm(list=ls())

library(splines)
library(survival)
library(MCMCpack)

n<-100
B<-1000

n.gran<-2
n.each<-50

vec<-c(4,5,6)


for(k in 1:3){
  for(v in 1:3){	
    
    E2<-rep(sort(c(1:7,22:36,8:21,8:21)),each=n.gran)
    
    E<-E2
    
    E_t<-rep(sort(c(1:7,22:36,8:21,8:21)),each=200)
    
    times.ev<-sort(unique(E))
    
    le.est<-c.est<-le.est.naive<-t.est<-matrix(0,nrow=length(times.ev),ncol=B)
    
    
    for(i in 1:B){
      
      
      
      if(k==1){
        # 16
        l.param<-log(E2,base=1.5)+1.3
        l.param.inv<-(sum(range(l.param))+0.0-l.param)-min(0,range(l.param)[1])+0.05
        
        #l.param.inv<-(sum(range(l.param))+0.0-l.param)
        
        l.param2<-(sum(range(l.param))+0-l.param.inv)-(mean(l.param)-mean(l.param.inv))
        
        l.param2<-0.60*(l.param2-min(0,range(l.param2)[1])+0.05)
        
        
        L<-rgamma(n.gran*n.each,l.param.inv,1/l.param.inv)
        
        
        L2_16_pos<-rgamma(n.gran*n.each,l.param2,1/l.param2)
        L_t<-rgamma(10000,l.param2,1/l.param2)
        
      }
      #L2<-rgamma(n.gran*n.each,l.param,1/l.param)
      
      
      if(k==2){
        #25
        l.param<-log(E2,base=1.5)+2.5
        l.param.inv<-(sum(range(l.param))+0-l.param)
        l.param2<-(sum(range(l.param))+0-l.param.inv)-(mean(l.param)-mean(l.param.inv))
        
        l.param2<-0.75*(l.param2-min(0,range(l.param2)[1])+0.05)
        
        L<-rgamma(n.gran*n.each,l.param.inv,1/l.param.inv)
        L2_25_pos<-rgamma(n.gran*n.each,l.param2,1/l.param2)
        L_t<-rgamma(10000,l.param2,1/l.param2)
        
      }
      
      
      if(k==3){
        #36
        l.param<-log(E2,base=1.5)+3.4
        l.param.inv<-(sum(range(l.param))+0-l.param)
        l.param2<-(sum(range(l.param))+0-l.param.inv)-(mean(l.param)-mean(l.param.inv))
        
        l.param2<-0.88*(l.param2-min(0,range(l.param2)[1])+0.05)
        
        L<-rgamma(n.gran*n.each,l.param.inv,1/l.param.inv)
        
        L2_36_pos<-rgamma(n.gran*n.each,l.param2,1/l.param2)
        L_t<-rgamma(10000,l.param2,1/l.param2)
        
      }
      X<-rgamma(n.gran*n.each,vec[v],1/vec[v])
      X_t = rgamma(10000,vec[v],1/vec[v])
      
      cen<-ifelse(E<=L,E,L)
      y<-ifelse(X<=cen,X,cen)
      
      e.cen<-as.numeric(E==y)
      c.cen<-as.numeric(cen==y)
      l.cen<-as.numeric(L==y)
      x.cen<-as.numeric(X==y)
      
      kmc<-survfit(Surv(y,c.cen)~1)
      kml<-survfit(Surv(y,l.cen)~E)
      kme<-survfit(Surv(y,e.cen)~1)
      kme.cdf<-survfit(Surv(E,rep(1,length(E)))~1)
      kmcen.cdf<-survfit(Surv(cen,rep(1,length(cen)))~1)
      
      sum.c<-summary(kmc,times=times.ev,extend=TRUE)
      sum.e<-summary(kme,times=times.ev,extend=TRUE)
      sum.e.cdf<-summary(kme.cdf,times=times.ev,extend=TRUE)
      sum.cen.cdf<-summary(kmcen.cdf,times=times.ev,extend=TRUE)
      
      
      l.surv<-matrix(0,nrow=length(times.ev),ncol=(length(times.ev)-1))
      
      for(j in 1:(length(times.ev)-1)){
        inds<-which(E>times.ev[j])
        kml.nest<-survfit(Surv(y[inds],l.cen[inds])~1)
        l.surv[,j]<-summary(kml.nest,times=times.ev,extend=TRUE)$surv
      }
      
      l.naive<-summary(survfit(Surv(y,l.cen)~1),times=times.ev,extend=TRUE)$surv
      
      
      c.est[,i]<-sum.c$surv
      le.est[,i]<-c(diag(l.surv),0)*sum.e.cdf$surv
      le.est.naive[,i]<-l.naive*sum.e.cdf$surv
      t.est[,i] <- sum.cen.cdf$surv
      
      
      #print(i)
      
    }
    
    
    
    
    library(ggplot2)
    library(plyr)
    
    
    le.curve<-rowMeans(le.est)
    c.curve<-rowMeans(c.est)
    t.curve<-rowMeans(t.est) 
    
    fun.var = function(x){
      var(x)
    }
    
    var.c.est<-aaply(c.est,1, fun.var)
    var.le.est<-aaply(le.est,1,fun.var)
    
    fun.q <- function(x){
      quantile(x, c(0.025, 0.975))
    }
    
    quant.c = aaply(c.est,1,fun.q)
    quant.le = aaply(le.est,1,fun.q)
    
    inds<-which(var.c.est>var.le.est)
    r.eff<-range(times.ev[inds])
    
    dist.dat<-data.frame(cbind(E = E_t,L = L_t, X = X_t))
    
    dat <- data.frame(
      xmin = r.eff[1],
      xmax = r.eff[2],
      ymin = 0,
      ymax = 0.20
    )
    
    thru<-40
    
    d.var<-data.frame(
      times.ev= times.ev[1:thru],
      le.var= var.le.est[1:thru],
      c.var= var.c.est[1:thru],
      quo= (var.le.est/var.c.est)[1:thru])
    
    
    d.surv<-data.frame(
      surv=le.curve,
      surv.c = c.curve,
      time=times.ev
    )
    
    d.surv.short<-data.frame(
      surv=le.curve[1:thru],
      surv.c = c.curve[1:thru],
      surv.t = t.curve[1:thru],
      time=times.ev[1:thru]
    )
    
    d.quant <- data.frame(
      times.ev = times.ev[1:thru],
      le.lower = as.numeric(quant.le[,1][1:thru]),
      le.upper = as.numeric(quant.le[,2][1:thru]),
      c.lower = as.numeric(quant.c[,1][1:thru]),
      c.upper = as.numeric(quant.c[,2][1:thru])
      
    )
    
    p1<-ggplot()+ coord_cartesian(x=c(0,90),y=c(0,0.09))+
      geom_density(data=dist.dat,aes(x=E),fill="#333333",alpha=0.7)+
      geom_density(data=dist.dat,aes(x=L),fill="#AAAAAA",alpha=0.7)+
      geom_density(data=dist.dat,aes(x=X),fill="#EEEEEE",alpha=0.7)+
      # geom_rect(data=dat,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),color="blue",fill="orange",alpha=0.3)+
      ylab("density")+xlab("time")+
      labs(title = "a)")+
      theme(plot.title = element_text(hjust = -0.12),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14))
    
  
    p3<-ggplot(data=data.frame())+geom_line(data=d.var,aes(x=times.ev,y=quo),colour="#555555",size=1.8,linetype=1)+
      coord_cartesian(x=c(0,70),y=c(0,1))+
      xlab("time")+ylab("rel. eff.")+
      geom_hline(yintercept=1,linetype=2)+
      labs(title = "b)")+
      theme(plot.title = element_text(hjust = -0.12),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14))
    # print(p3)
    
    p4<-ggplot(data=data.frame())+
      coord_cartesian(x=c(0,70),y=c(0,1))+xlab("time")+
      ylab("survival")+
      geom_line(data=d.surv.short,aes(x=time,y=surv.c),colour="#999999",size=.5,linetype=1)+
      geom_line(data=d.surv.short,aes(x=time,y=surv),colour=1,size=.5,linetype=3)+
      geom_line(data=d.surv.short,aes(x=time,y=surv.t),colour="red",size=.5,linetype=2)+
      geom_ribbon(data = d.quant, aes(x = times.ev, ymin = le.lower, ymax = le.upper), alpha =0.6)+
      geom_ribbon(data = d.quant, aes(x = times.ev, ymin = c.lower, ymax = c.upper), alpha =0.3)+
      labs(title = "c)")+
      theme(plot.title = element_text(hjust = -0.12),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14))
    
    
    
    
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      require(grid)
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
        print(plots[[1]])
        
      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
          # Get the i,j matrix positions of the regions that contain this subplot
          matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
          
          print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                          layout.pos.col = matchidx$col))
        }
      }
    }
    
    
    b<-gsub("q",(100*v+10*k+0),"graph_neg_q.pdf")
    pdf(b,width=5,heigh=7)
    multiplot(p1,p3,p4,cols=1)
    dev.off()
    
    print(v)
    print(k)
    
  }	
}
