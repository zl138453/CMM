CMM=function(data,outcome,med,pred,cov_con=NULL,cov_cat=NULL,weight=NULL,family="identity",boot=5000){

  if(family=="identity"){
    result=md_linear(data,outcome,med,pred,cov_con=cov_con,cov_cat=cov_cat,weight=weight,boot=boot)
  }else{
    result=md_logistic(data,outcome,med,pred,cov_con=cov_con,cov_cat=cov_cat,weight=weight,boot=boot)
  }

}

###m is coefficient vector
eest1=function(m){
  m=c(m,0)
  rsa=matrix(NA,1,length(m))

  for(i in 1:(length(m))){
    if(i==1){
      rsa[1,i]=m[i]*sqrt((length(m)-i)/(length(m)+1-i))
    }else{
      b=0
      for(j in 1:(i-1)){
        b=b+m[j]*sqrt(1/((length(m)+1-j)*(length(m)-j)))
      }
      rsa[1,i]=m[i]*sqrt((length(m)-i)/(length(m)+1-i))-b

    }


  }

  return(rsa)
}


md_linear=function(data,outcome,med,pred,cov_con=NULL,cov_cat=NULL,weight=NULL,boot=5000){
  nmed=length(med)
  n1=ncol(data)

  ###Deal with categorical

  if(is.null(cov_cat)){
    cov=cov_con
  }else{

    for( i in cov_cat){
      na=ncol(data)
      name1=colnames(data)[i]
      name2=names(table(data[,i]))
      data=cbind(data,dummy_cols(data[,i])[,-c(1,2)])
      nb=ncol(data)
      colnames(data)[(na+1):nb]=paste(name1,"_",name2[2:length(name2)],sep="")

    }
    n2=ncol(data)
    ###All coveriate variables
    cov=c((n1+1):n2,cov_con)

  }



  ###ILR transform

  df  <- data[,med]


  df1 <- pivotCoord(df,pivotvar=1)

  colnames(df1) <- gsub("-", ".", colnames(df1))   # swap - for . in variable names

  data <- data.frame(data,df1)

  n3=ncol(data)

  ########Create matrix

  ###ilr matrix
  ilr=as.matrix(data[,(n3-nmed+2):n3])

  ###Non ilr matrix
  com=as.matrix(data[,med])

  ###Covariate
  if(is.null(cov)){
    cov1=NULL

  }else{
    cov1=as.matrix(data[,cov])
  }


  ###Predictor
  pred1=as.matrix(data[,pred])

  ####Outcome
  out=as.matrix(data[,outcome])

  ###weight

  #w=ifelse(is.null(weight),as.matrix(rep(1,nrow(data))),as.matrix(data[,weight]))

  if(is.null(weight)){
    w=as.matrix(rep(1,nrow(data)))
  }else{
    w=as.matrix(data[,weight])
  }


  est1=matrix(NA,6,nmed)

  for(i in 1:nmed){
    if(is.null(cov1)){
      a=summary(lm(log(com[,i])~pred1,weights=w))$coefficients

    }else{
      a=summary(lm(log(com[,i])~pred1+cov1,weights=w))$coefficients
    }

    est1[1,i]=a[2,1]
    est1[2,i]=a[2,4]
  }
  colnames(est1)=colnames(data[,med])

  est2=matrix(NA,4,nmed)

  if(is.null(cov1)){
    b=summary(lm(out~pred1+ilr,weights=w))$coefficients

  }else{
    b=summary(lm(out~pred1+ilr+cov1,weights=w))$coefficients
  }

  est2[1,1:(nmed-1)]=b[3:(1+nmed),1]
  est2[2,1:(nmed-1)]=b[3:(1+nmed),4]
  est2[1,nmed]=b[2,1]
  est2[2,nmed]=b[2,4]

  est3=eest1(est2[1,1:(nmed-1)])

  ###Need to recode to automatically rename

  colnames(est3)=colnames(data[,med])

  ef=matrix(NA,5,nmed+3)
  ef[1,1:nmed]=est1[1,]*est3

  ##total indirect ef
  ef[1,1+nmed]=sum(ef[1,1:nmed])
  ##Total effect
  ef[1,2+nmed]=est2[1,nmed]
  ef[1,3+nmed]=ef[1,1+nmed]+ef[1,2+nmed]


  ###bootstrp
  est1b=matrix(NA,boot,nmed)
  colnames(est1b)=colnames(data[,med])
  est2b=matrix(NA,boot,nmed)
  efb=matrix(NA,boot,nmed+3)

  for(j in 1:boot){
    s=sample(1:nrow(data),nrow(data),replace = T)

    ilr1=ilr[s,]

    ###Non ilr matrix
    com1=com[s,]

    ###Covariate
    cov2=cov1[s,]

    ###Predictor
    pred2=pred1[s]

    ####Outcome
    out1=out[s]

    ###weight
    w1=w[s]


    for(i in 1:nmed){
      if(is.null(cov2)){
        a=summary(lm(log(com1[,i])~pred2,weights=w1))$coefficients
      }else{
        a=summary(lm(log(com1[,i])~pred2+cov2,weights=w1))$coefficients
      }

      est1b[j,i]=a[2,1]
    }


    if(is.null(cov2)){
      b=summary(lm(out1~pred2+ilr1,weights=w1))$coefficients
    }else{
      b=summary(lm(out1~pred2+ilr1+cov2,weights=w1))$coefficients

    }

    est2b[j,1:(nmed-1)]=b[3:(1+nmed),1]
    est2b[j,nmed]=b[2,1]

    est3=eest1(est2b[j,1:(nmed-1)])


    colnames(est3)=colnames(data[,med])

    efb[j,1:nmed]=est1b[j,]*est3

    ##total indirect ef
    efb[j,1+nmed]=sum(efb[j,1:nmed])
    ##Direct effect
    efb[j,2+nmed]=est2b[j,nmed]
    ##Total effect
    efb[j,3+nmed]=efb[j,1+nmed]+ efb[j,2+nmed]



  }

  m1=apply(est1b,2,mean)
  s1=apply(est1b,2,sd)
  lp1=apply(est1b,2,quantile,prob=0.025)
  up1=apply(est1b,2,quantile,prob=0.975)
  est1[3,]=m1
  est1[4,]=s1
  est1[5,]=lp1
  est1[6,]=up1

  m2=apply(efb,2,mean)
  s2=apply(efb,2,sd)
  lp2=apply(efb,2,quantile,prob=0.025)
  up2=apply(efb,2,quantile,prob=0.975)
  ef[2,]=m2
  ef[3,]=s2
  ef[4,]=lp2
  ef[5,]=up2

  rownames(est1)=c("Est","Est_SD","boot_Mean","boot_SD","boot_lb","boot_ub")
  rownames(ef)=c("Est","boot_Mean","boot_SD","boot_lb","boot_ub")
  colnames(ef)=c(paste(colnames(data[,med]),"_IE",sep=""),"TIE","DE","TE")

  ###Contribution to total effect
  ppa=efb
  ppa[,1:nmed]=ppa[,1:nmed]/ppa[,3+nmed]
  ppa[,1+nmed]=ppa[,1+nmed]/ppa[,3+nmed]
  ppa[,2+nmed]=ppa[,2+nmed]/ppa[,3+nmed]
  ppa=ppa[,1:(2+nmed)]

  pp=matrix(NA,4,2+nmed)
  pp[1,]=apply(ppa,2,mean)
  pp[2,]=apply(ppa,2,sd)
  pp[3,]=apply(ppa,2,quantile,prob=0.025)
  pp[4,]=apply(ppa,2,quantile,prob=0.975)

  colnames(pp)=c(paste(colnames(data[,med]),"_IE",sep=""),"TIE","DE")
  rownames(pp)=c("boot_Mean","boot_SD","boot_lb","boot_ub")

  gg<-data.frame(cbind(name=colnames(pp),t(pp)))
  gg[,2:ncol(gg)]<-sapply( gg[,2:ncol(gg)], as.numeric)

  fig<-gg %>%
    mutate(name = fct_reorder(name, desc(1:nrow(gg)))) %>%
    ggplot( aes(x=name, y=boot_Mean)) +
    geom_pointrange(aes(ymin=boot_lb, ymax=boot_ub))+theme_bw()+xlab("")+ylab("Relative Effects")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_y_continuous(breaks = seq(floor(min(gg$boot_lb)*10)/10, ceiling(max(gg$boot_ub)*10)/10,
                                    by = ceiling((ceiling(max(gg$boot_ub)*10)/10-floor(min(gg$boot_lb)*10)/10))/10))+
    coord_flip()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  IDE=ef[c(2,4,5),-c((ncol(ef)-1):ncol(ef))]
  DE=ef[c(2,4,5),ncol(ef)-1]
  TE=ef[c(2,4,5),ncol(ef)]


  gg1<-data.frame(cbind(name=colnames(ef),t(ef)))
  gg1[,2:ncol(gg1)]<-sapply(gg1[,2:ncol(gg1)], as.numeric)

  fig1<-gg1 %>%
    mutate(name = fct_reorder(name, desc(1:nrow(gg1)))) %>%
    ggplot( aes(x=name, y=boot_Mean)) +
    geom_pointrange(aes(ymin=boot_lb, ymax=boot_ub))+theme_bw()+xlab("")+ylab("Mediation Effect (Coefficient)")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    coord_flip()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))



  mylist=list("Indirect.effect"=IDE,"Direct.effect"=DE,"Total.effect"=TE,"Mediation.effect.plot"=fig1,"Relative.Effects.plot"=fig)
  mylist

}

md_logistic=function(data,outcome,med,pred,cov_con=NULL,cov_cat=NULL,weight=NULL,boot=5000){
  nmed=length(med)
  n1=ncol(data)

  ###Deal with categorical

  if(is.null(cov_cat)){
    cov=cov_con
  }else{

    for( i in cov_cat){
      na=ncol(data)
      name1=colnames(data)[i]
      name2=names(table(data[,i]))
      data=cbind(data,dummy_cols(data[,i])[,-c(1,2)])
      nb=ncol(data)
      colnames(data)[(na+1):nb]=paste(name1,"_",name2[2:length(name2)],sep="")

    }
    n2=ncol(data)
    ###All coveriate variables
    cov=c((n1+1):n2,cov_con)

  }



  ###ILR transform

  df  <- data[,med]


  df1 <- pivotCoord(df,pivotvar=1)

  colnames(df1) <- gsub("-", ".", colnames(df1))   # swap - for . in variable names

  data <- data.frame(data,df1)

  n3=ncol(data)

  ########Create matrix

  ###ilr matrix
  ilr=as.matrix(data[,(n3-nmed+2):n3])

  ###Non ilr matrix
  com=as.matrix(data[,med])

  ###Covariate
  if(is.null(cov)){
    cov1=NULL

  }else{
    cov1=as.matrix(data[,cov])
  }


  ###Predictor
  pred1=as.matrix(data[,pred])

  ####Outcome
  out=as.matrix(data[,outcome])

  ###weight

  #w=ifelse(is.null(weight),as.matrix(rep(1,nrow(data))),as.matrix(data[,weight]))

  if(is.null(weight)){
    w=as.matrix(rep(1,nrow(data)))
  }else{
    w=as.matrix(data[,weight])
  }


  est1=matrix(NA,6,nmed)

  for(i in 1:nmed){
    if(is.null(cov1)){
      a=summary(lm(log(com[,i])~pred1,weights=w))$coefficients

    }else{
      a=summary(lm(log(com[,i])~pred1+cov1,weights=w))$coefficients
    }

    est1[1,i]=a[2,1]
    est1[2,i]=a[2,4]
  }
  colnames(est1)=colnames(data[,med])

  est2=matrix(NA,4,nmed)

  if(is.null(cov1)){
    uuu=svydesign(ids=~1, weights=~w,data=data.frame(w=w, out=out,pred1=pred1,ilr))
    kk=as.formula(paste(colnames(uuu$variables)[2],paste(colnames(uuu$variables)[3:ncol(uuu$variables)],collapse ="+"),sep="~"))
    b=summary(svyglm(kk,design=uuu, family=binomial))$coefficients
  }else{
    uuu=svydesign(ids=~1, weights=~w,data=data.frame(w=w, out=out,pred1=pred1,ilr,cov1))
    kk=as.formula(paste(colnames(uuu$variables)[2],paste(colnames(uuu$variables)[3:ncol(uuu$variables)],collapse ="+"),sep="~"))
    b=summary(svyglm(kk,design=uuu, family=binomial))$coefficients
  }

  est2[1,1:(nmed-1)]=b[3:(1+nmed),1]
  est2[2,1:(nmed-1)]=b[3:(1+nmed),4]
  est2[1,nmed]=b[2,1]
  est2[2,nmed]=b[2,4]

  est3=eest1(est2[1,1:(nmed-1)])

  ###Need to recode to automatically rename

  colnames(est3)=colnames(data[,med])

  ef=matrix(NA,5,nmed+3)
  ef[1,1:nmed]=exp(est1[1,]*est3)

  ##total indirect ef
  ef[1,1+nmed]=prod(ef[1,1:nmed])
  ##Total effect
  ef[1,2+nmed]=exp(est2[1,nmed])
  ef[1,3+nmed]=ef[1,1+nmed]*ef[1,2+nmed]


  ###bootstrp
  est1b=matrix(NA,boot,nmed)
  colnames(est1b)=colnames(data[,med])
  est2b=matrix(NA,boot,nmed)
  efb=matrix(NA,boot,nmed+3)

  for(j in 1:boot){
    s=sample(1:nrow(data),nrow(data),replace = T)

    ilr1=ilr[s,]

    ###Non ilr matrix
    com1=com[s,]

    ###Covariate
    cov2=cov1[s,]

    ###Predictor
    pred2=pred1[s]

    ####Outcome
    out1=out[s]

    ###weight
    w1=w[s]


    for(i in 1:nmed){
      if(is.null(cov2)){
        a=summary(lm(log(com1[,i])~pred2,weights=w1))$coefficients
      }else{
        a=summary(lm(log(com1[,i])~pred2+cov2,weights=w1))$coefficients
      }

      est1b[j,i]=a[2,1]
    }


    if(is.null(cov2)){
      uuu=svydesign(ids=~1, weights=~w1,data=data.frame(w1=w1, out1=out1,pred2=pred2,ilr1))
      kk=as.formula(paste(colnames(uuu$variables)[2],paste(colnames(uuu$variables)[3:ncol(uuu$variables)],collapse ="+"),sep="~"))
      b=summary(svyglm(kk,design=uuu, family=binomial))$coefficients

    }else{
      uuu=svydesign(ids=~1, weights=~w1,data=data.frame(w1=w1, out1=out1,pred2=pred2,ilr1,cov2))
      kk=as.formula(paste(colnames(uuu$variables)[2],paste(colnames(uuu$variables)[3:ncol(uuu$variables)],collapse ="+"),sep="~"))
      b=summary(svyglm(kk,design=uuu, family=binomial))$coefficients

    }

    est2b[j,1:(nmed-1)]=b[3:(1+nmed),1]
    est2b[j,nmed]=b[2,1]

    est3=eest1(est2b[j,1:(nmed-1)])


    colnames(est3)=colnames(data[,med])

    efb[j,1:nmed]=exp(est1b[j,]*est3)

    ##total indirect ef
    efb[j,1+nmed]=prod(efb[j,1:nmed])
    ##Direct effect
    efb[j,2+nmed]=exp(est2b[j,nmed])
    ##Total effect
    efb[j,3+nmed]=efb[j,1+nmed]*efb[j,2+nmed]



  }

  m1=apply(est1b,2,mean)
  s1=apply(est1b,2,sd)
  lp1=apply(est1b,2,quantile,prob=0.025)
  up1=apply(est1b,2,quantile,prob=0.975)
  est1[3,]=m1
  est1[4,]=s1
  est1[5,]=lp1
  est1[6,]=up1

  m2=exp(apply(log(efb),2,mean))
  s2=exp(apply(log(efb),2,sd))
  lp2=exp(apply(log(efb),2,quantile,prob=0.025))
  up2=exp(apply(log(efb),2,quantile,prob=0.975))
  ef[2,]=m2
  ef[3,]=s2
  ef[4,]=lp2
  ef[5,]=up2

  rownames(est1)=c("Est","Est_SD","boot_Mean","boot_SD","boot_lb","boot_ub")
  rownames(ef)=c("Est","boot_Mean","boot_SD","boot_lb","boot_ub")
  colnames(ef)=c(paste(colnames(data[,med]),"_IE",sep=""),"TIE","DE","TE")

  ###Contribution to total effect
  ppa=efb
  ppa[,1:nmed]=log(ppa[,1:nmed])/log(ppa[,3+nmed])
  ppa[,1+nmed]=log(ppa[,1+nmed])/log(ppa[,3+nmed])
  ppa[,2+nmed]=log(ppa[,2+nmed])/log(ppa[,3+nmed])
  ppa=ppa[,1:(2+nmed)]

  pp=matrix(NA,4,2+nmed)
  pp[1,]=apply(ppa,2,mean)
  pp[2,]=apply(ppa,2,sd)
  pp[3,]=apply(ppa,2,quantile,prob=0.025)
  pp[4,]=apply(ppa,2,quantile,prob=0.975)

  colnames(pp)=c(paste(colnames(data[,med]),"_IE",sep=""),"TIE","DE")
  rownames(pp)=c("boot_Mean","boot_SD","boot_lb","boot_ub")

  gg<-data.frame(cbind(name=colnames(pp),t(pp)))
  gg[,2:ncol(gg)]<-sapply( gg[,2:ncol(gg)], as.numeric)

  fig<-gg %>%
    mutate(name = fct_reorder(name, desc(1:nrow(gg)))) %>%
    ggplot( aes(x=name, y=boot_Mean)) +
    geom_pointrange(aes(ymin=boot_lb, ymax=boot_ub))+theme_bw()+xlab("")+ylab("Relative Effects")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_y_continuous(breaks = seq(floor(min(gg$boot_lb)*10)/10, ceiling(max(gg$boot_ub)*10)/10,
                                    by = ceiling((ceiling(max(gg$boot_ub)*10)/10-floor(min(gg$boot_lb)*10)/10))/10))+
    coord_flip()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  IDE=ef[c(2,4,5),-c((ncol(ef)-1):ncol(ef))]
  DE=ef[c(2,4,5),ncol(ef)-1]
  TE=ef[c(2,4,5),ncol(ef)]


  gg1<-data.frame(cbind(name=colnames(ef),t(ef)))
  gg1[,2:ncol(gg1)]<-sapply(gg1[,2:ncol(gg1)], as.numeric)

  fig1<-gg1 %>%
    mutate(name = fct_reorder(name, desc(1:nrow(gg1)))) %>%
    ggplot( aes(x=name, y=boot_Mean)) +
    geom_pointrange(aes(ymin=boot_lb, ymax=boot_ub))+theme_bw()+xlab("")+ylab("Mediation Effect (OR)")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    coord_flip()+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))



  mylist=list("Indirect.effect"=IDE,"Direct.effect"=DE,"Total.effect"=TE,"Mediation.effect.plot"=fig1,"Relative.Effects.plot"=fig)
  mylist

}





