## A function for selecting variables that are important predictors in a categorical
  # (can, but does not have to be binary) random forest model.

## Leverages the conditional variable importance from 'cforest' in the 'party' package
  # to extract more robust variable importance rankings.

## For binary variables, calculates conditional mean decrease accuracy on AUC measure
  # rather than % accuracy.

## In order to save time/computing resources that can be prohibitive for cforest, data
  # are split, and sent to cluster nodes.

## Variable importance metrics from each cforest from each cluster node are averaged.

## To avoid splitting data into multiple nodes, set ncpus to 1.

## Variable selection proceeds with reverse selection, in a proportional stepwise function
  # according to the setting of 'ndrop' (see below).

## The tuning parameter, 'keep.prop' determines how many variables are kept after the selection procedure has
  # created a wide array of models to compare.

## The tuning parameter 'mb' sets the minimum size for terminal nodes in the cforest models.  Keeping this small
  # tends to lead to models that produce stronger signals in terms of ranking variable importance, but comes
  # with a processing-time cost.
  # It is worthwhile to test this function with larger values (> or equal to the number of observations in the smallest category).

## 'stop.at' is a tuning parameter that can save time.  It stops the stepwise procedure when the current model accuracy drops below the maximum model accuracy by this threshold.

## The series of boxplots produced show variable importance from the sub-models, and the accuracy of each model at each step.

## When stepwise variable selection has completed, the mean accuracy statistics for each model at each step are plotted, and a lowess curve fitted.
  # The final list of variables selected correspond with those included in the models driving this last plot.

## package dependencies: snowfall, party and caTools must be installed.

VarSel_cforest<-function(dat, ## Data frame containing y-variable as the first column, x-variables to reduce follow.
                         ndrop = 2, ## number of variables to drop at each step, to begin with. The number of variables to drop at each step becomes smaller as there are fewer variables to consider.  If this set to less than 1,  proportion of variables to crop at each step. If this is set to greater than 1, a proportion is calculated according to the number of variables selected, and that proportion is used to determine how many variables are dropped at each step.
                         keep.prop = .9, ## Cutoff determining which variables to keep in the final selection. All variables included in the models with an accuracy greater than keep.prop times the maximum accuracy are returned by the function. This thresholding decision comes after stepwise variable selection has completed.
                         mb = 20, ## minimum size for terminal nodes in cforest models used to get variable importance.
                         stop.at = .005, ## Stop variable selection loop when maximum accuracy is greater than current accuracy by this number.
                         ncpus = 6 ## number of cpus to use in the cluster.
                         ){
  require(snowfall)
  sfInit(T,ncpus)
  sfLibrary(party)
  cv<-colnames(dat)[2:ncol(dat)]
  sfExport(list =c("mb","cv"))
  vsl<-list()
  i<-0
  if(ndrop >=1){nd<-ndrop;ndrop<-ndrop/length(cv)
  }else{nd<-max(c(1,round(ndrop * length(cv))))}
  cat(paste("one\nndrop:",ndrop,"nd:",nd ,"\n"))
  minacc<-0
  while(length(cv)>1 & minacc < stop.at){
    cat("\n",length(cv),"\n")
    vsl[[i<-i+1]]<-VarSel_cforest_GetOne(x=dat[,c(1,which(colnames(dat) %in% cv))],samp = NULL, mb = mb, ncpus = ncpus)
    cv<-cv[!cv %in%  names(vsl[[i]][[1]])[1:nd]]
    nd<-round(ndrop*(length(cv)))
    if(nd < 1)nd<-1
    cat(paste("ndrop:",ndrop,"nd:",nd,"\n"))
    if(length(vsl) > 5){
      accs<-do.call(rbind,lapply(vsl,function(x){return(x[[2]])}))
      l1<-lowess(accs,f = .8)
      lastacc<-l1$y[length(l1$y)]
      maxacc<-max(l1$y)
      minacc<-maxacc - lastacc
      print(cv)
      cat("minacc:",minacc,"\nlastacc:",lastacc,"\nmaxacc:",maxacc,"\nstop.at:",stop.at,"\n")
    }
  }
  sfStop()
  accs<-do.call(rbind,lapply(vsl,function(x){return(x[[2]])}))
  par(mfrow =c(ncol(accs),1))
  plot(accs)
  if(length(accs)< 6){lines(l1<-lowess(accs,f = .8),lty = 2)
  }else{lines(l1<-lowess(accs),lty = 2)}
  browser()
  co<-max(l1$y)*keep.prop
  wm<-max(1,max(which (l1$y > co))-1)## selects final model as the one before the step identified by the keep.prop parameter.
                                     ## This is often the second-to-last model.

  vl<-names(vsl[[wm]][[1]])
  return(vl)
}

## This is an internal function called by VarSel_cforest
VarSel_cforest_GetOne<-function(x,
                                samp = NULL,
                                mb = nrow(x)/100,
                                ncpus = 6)## Note: for this last parameter, it would be good to recode this to simply match the number of cluster nodes created in the primary function's call to sfInit.
{
  if(!is.null(samp))x<-x[sample(1:nrow(x),samp),]
  ind<-rep(1:ncpus,length.out = nrow(x))
  x<-split(x,ind)
  vil<-sfLapply(x,function(z2){
    if(lt1<-length(table(z2$yv<-factor(z2$yv)))>1){
      cf<-c()
      counter<-1
      while(!"RandomForest" %in% class(cf) & counter < 4){
        cf<-try({party::cforest(yv~.,data = z2,
                                controls =cforest_unbiased(minbucket = mb))
        })
        counter<-counter+1
      }
      if("RandomForest" %in% class(cf)){
        if(lt1 == 2){
          ## calculate varimp, and accuracy based on AUC
          print(table(z2$yv))
          pred1<-do.call(rbind,predict(cf,type = "prob"))
          acc<-caTools::colAUC(X = pred1, y= z2$yv, plotROC=F, alg=c("Wilcoxon","ROC"))
          vi<-party::varimpAUC(cf,conditional = T, OOB = T)
          return(list(vi,acc[2]))
        }else{
          ## " calculate ordinary varimp, accuracy based on AUC"
          pred1<-do.call(rbind,predict(cf,type = "prob"))
          acc<-mean(caTools::colAUC(X = pred1, y= z2$yv, plotROC=F, alg=c("Wilcoxon","ROC")))
          vi<-party::varimp(cf,conditional = T, OOB = T) ## Note: for multi-category RF models, varimp is based on % accuracy, rather than AUC as it is for binary models.
          return(list(vi,acc))
        }
      }else{return("cforest failed")}
    }else{return("only one category in yvar")}
  })
  vil1<-do.call(cbind,lapply(vil,function(x){if(class(x) == "character"){return(NULL)}else{return(x[[1]])}}))
  vil2<-do.call(rbind,lapply(vil,function(x){if(class(x) == "character"){return(NULL)}else{return(x[[2]])}}))
  boxplot(t(vil1[order(vil3<-apply(vil1,1,mean)),]),las = 2)
  acc<-mean(vil2)
  mtext(paste("   ",round(acc,3)),adj = 0,side = 3, line = -3)
  return(list(sort(vil3),acc))
}
