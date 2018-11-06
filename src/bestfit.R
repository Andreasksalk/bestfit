### The function that does all the calculations, it can take some seconds to run
bestfit <- function(datcolumn,bins){
  
  fitgood = data.table(distribution = c("normal","exponential","poisson","lognormal","weibull","gamma","triangular","uniform"),sqerror = c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0))
  x = datcolumn
  xmin = (min(x))
  xmax = (max(x))
  
  
  res <- hist(x,plot=FALSE, breaks = seq(xmin,xmax,l=bins+1))
  divM = sum(res$density)
  
  # Square error of normal distribution
  ft.norm <- fitdistr(x,"normal")
  x.norm = vector()
  for(i in 1:(length(res$breaks)-1)){
    x.norm[i] = pnorm(res$breaks[i+1],ft.norm$estimate[1],ft.norm$estimate[2])-pnorm(res$breaks[i],ft.norm$estimate[1],ft.norm$estimate[2])
  }
  
  chi.norm = chisq.test(p = res$density/divM, x = x.norm*length(x))
  ks.norm = ks.test(x = x , y = rnorm(x,ft.norm$estimate[1],ft.norm$estimate[2]))
  sqer.norm = sum((res$density/divM-x.norm)^2)
  sum(((res$density/divM-x.norm)^2))
  
  fitgood[1,2] = sqer.norm
  
  #Square error of exponential distribution
  ft.expo <- fitdistr(x,"exponential")
  x.expo = vector()
  for(i in 1:(length(res$breaks)-1)){
    x.expo[i]<-pexp(res$breaks[i+1], ft.expo$estimate[1])-pexp(res$breaks[i],ft.expo$estimate[1])
  }
  length(res$breaks)
  
  
  chi.exp = chisq.test(p = res$density/divM,x = x.expo*length(x))
  
  ks.exp = ks.test(x, rexp(x,ft.expo$estimate[1]))
  sqer.exp = sum((res$density/divM-x.expo)^2)
  fitgood[2,2] = sqer.exp
  
  #Square error of poisson (does not makes sense on continously data)
  ft.pois <- fitdistr(x,"poisson")
  
  x.pois= vector()
  for(i in 1:(length(res$breaks)-1)){
    x.pois[i] = ppois(res$breaks[i+1], ft.pois$estimate[1])-ppois(res$breaks[i],ft.pois$estimate[1])
  }
  
  chi.pois = chisq.test(p = res$density/divM,x = x.pois*length(x))
  ks.pois = ks.test(x, rpois(x,ft.pois$estimate[1]))
  sqer.pois = sum((res$density/divM-x.pois)^2)
  fitgood[3,2] = sqer.pois
  
  
  #Square error of Lognormal
  ft.lognorm <- fitdistr(x,"lognormal")
  
  x.lognorm = vector()
  for(i in 1:(length(res$breaks)-1)){
    x.lognorm[i] <-plnorm(res$breaks[i+1],ft.lognorm$estimate[1],ft.lognorm$estimate[2])-plnorm(res$breaks[i],ft.lognorm$estimate[1],ft.lognorm$estimate[2])
  }
  
  chi.lognorm = chisq.test(p = res$density/divM,x = x.lognorm*length(x))
  ks.lognorm = ks.test(x, rlnorm(x,ft.lognorm$estimate[1],ft.lognorm$estimate[2]))
  sqer.lognorm = sum((res$density/divM-x.lognorm)^2)
  fitgood[4,2] = sqer.lognorm
  
  #Square error of Weibull
  ft.wei <- fitdistr(x,"weibull")
  
  x.wei = vector ()
  for(i in 1:(length(res$breaks)-1)){
    x.wei[i] <- pweibull(res$breaks[i+1], ft.wei$estimate[1],ft.wei$estimate[2])-pweibull(res$breaks[i], ft.wei$estimate[1],ft.wei$estimate[2])
  }
  
  chi.wei = chisq.test(p = res$density/divM,x = x.wei*length(x))
  ks.wei = ks.test(x, rweibull(x,ft.wei$estimate[1],ft.wei$estimate[2]))
  sqer.wei = sum((res$density/divM-x.wei)^2)
  fitgood[5,2] = sqer.wei
  sum(x.wei)
  
  #Square error of Gamma
  ft.gamma <- fitdistr(x,"gamma")
  
  x.gamma = vector()
  for(i in 1:(length(res$breaks)-1)){
    x.gamma[i] <- pgamma(res$breaks[i+1],ft.gamma$estimate[1],ft.gamma$estimate[2])-pgamma(res$breaks[i],ft.gamma$estimate[1],ft.gamma$estimate[2])
  }
  
  res$density
  
  chi.gamma = chisq.test(p = res$density/divM,x = x.gamma*length(x))
  ks.gamma = ks.test(x, rgamma(x,ft.gamma$estimate[1],ft.gamma$estimate[2]))
  sqer.gamma = sum((res$density/divM-x.gamma)^2)
  fitgood[6,2] = sqer.gamma
  
  #Square error of Triangular
  
  x.triang = vector()
  for(i in 1:(length(res$breaks)-1)){
    x.triang[i] <- ptriangle(res$breaks[i+1],min(x),max(x),res$mids[which(res$counts == max(res$counts))])-ptriangle(res$breaks[i],min(x),max(x),res$mids[which(res$counts == max(res$counts))])
  }
  
  chi.triang = chisq.test(p = res$density/divM,x = x.triang*length(x))
  ks.triang = ks.test(x, rtriangle(x,min(x),max(x),res$mids[which(res$counts == max(res$counts))]))
  sqer.triang = sum((res$density/divM-x.triang)^2)
  fitgood[7,2] = sqer.triang
  
  #Square error of Uniform
  x.uni = vector()
  for(i in 1:(length(res$breaks)-1)){
    x.uni[i] <- punif(res$breaks[i+1],min(x),max(x))-punif(res$breaks[i],min(x),max(x))
  }
  
  chi.uni = chisq.test(p = res$density/divM,x = x.uni*length(x))
  ks.uni = ks.test(x, runif(x,min(x),max(x)))
  sqer.uni = sum((res$density/divM-x.uni)^2)
  fitgood[8,2] = sqer.uni
  
  fitgood1 = fitgood[order(sqerror)]
  
  #######################
  # Getting the chi-squared and ks-test tests into the dashboard
  xparam11 = "NO VALUE"
  xparam12 = "NO VALUE"
  xparam13 = "NO VALUE"
  xparam21 = "NO VALUE"
  xparam22 = "NO VALUE"
  xparam23 = "NO VALUE"
  xparam31 = "NO VALUE"
  xparam32 = "NO VALUE"
  xparam33 = "NO VALUE"
  
  if(fitgood1[1,1] == "normal"){
    xsq1 = chi.norm$statistic
    df1 = chi.norm$parameter
    p1 = chi.norm$p.value
    d1 = ks.norm$statistic
    ksp1 = ks.norm$p.value
    xparam11 = paste("Mean", ft.norm$estimate[1], sep = ": ")
    xparam12 = paste("SD", ft.norm$estimate[2], sep = ": ")
    
  } else if (fitgood1[1,1] == "triangular"){
    xsq1 = chi.triang$statistic
    df1 = chi.triang$parameter
    p1 = chi.triang$p.value
    d1 = ks.triang$statistic
    ksp1 = ks.triang$p.value
    xparam11 = paste("Min", min(x), sep = ": ")
    xparam12 = paste("Max", max(x), sep = ": ")
    xparam13 = paste("Mode", res$mids[which(res$counts == max(res$counts))], sep = ": ")
  } else if (fitgood1[1,1] == "exponential"){
    xsq1 = chi.exp$statistic
    df1 = chi.exp$parameter
    p1 = chi.exp$p.value
    d1 = ks.exp$statistic
    ksp1 = ks.exp$p.value
    xparam11 = paste("Rate", ft.expo$estimate[1], sep = ": ")
  } else if (fitgood1[1,1] == "poisson"){
    xsq1 = chi.pois$statistic
    df1 = chi.pois$parameter
    p1 = chi.pois$p.value
    d1 = ks.pois$statistic
    ksp1 = ks.pois$p.value
    xparam11 = paste("Rate", ft.pois$estimate[1], sep = ": ")
  } else if (fitgood1[1,1] == "lognormal"){
    xsq1 = chi.lognorm$statistic
    df1 = chi.lognorm$parameter
    p1 = chi.lognorm$p.value
    d1 = ks.lognorm$statistic
    ksp1 = ks.lognorm$p.value
    xparam11 = paste("Mean", ft.lognorm$estimate[1], sep = ": ")
    xparam12 = paste("SD", ft.lognorm$estimate[2], sep = ": ")
  } else if (fitgood1[1,1] == "weibull"){
    xsq1 = chi.wei$statistic
    df1 = chi.wei$parameter
    p1 = chi.wei$p.value
    d1 = ks.wei$statistic
    ksp1 = ks.wei$p.value
    xparam11 = paste("Scale", ft.wei$estimate[1], sep = ": ")
    xparam12 = paste("Shape", ft.wei$estimate[2], sep = ": ")
  } else if (fitgood1[1,1] == "gamma"){
    xsq1 = chi.gamma$statistic
    df1 = chi.gamma$parameter
    p1 = chi.gamma$p.value
    d1 = ks.gamma$statistic
    ksp1 = ks.gamma$p.value
    xparam11 = paste("Shape",ft.gamma$estimate[1], sep = ": ")
    xparam12 = paste("Scale",ft.gamma$estimate[2], sep = ": ")
  } else if (fitgood1[1,1] == "uniform"){
    xsq1 = chi.uni$statistic
    df1 = chi.uni$parameter
    p1 = chi.uni$p.value
    d1 = ks.uni$statistic
    ksp1 = ks.uni$p.value
    xparam11 = paste("Min", min(x), sep = ": ")
    xparam12 = paste("Max", max(x), sep = ": ")
  }
  
  if(fitgood1[2,1] == "normal"){
    xsq2 = chi.norm$statistic
    df2 = chi.norm$parameter
    p2 = chi.norm$p.value
    d2 = ks.norm$statistic
    ksp2 = ks.norm$p.value
    xparam21 = paste("Mean", ft.norm$estimate[1], sep = ": ")
    xparam22 = paste("SD", ft.norm$estimate[2], sep = ": ")
  } else if (fitgood1[2,1] == "triangular"){
    xsq2 = chi.triang$statistic
    df2 = chi.triang$parameter
    p2 = chi.triang$p.value
    d2 = ks.triang$statistic
    ksp2 = ks.triang$p.value
    xparam21 = paste("Min", min(x), sep = ": ")
    xparam22 = paste("Max", max(x), sep = ": ")
    xparam23 = paste("Mode", res$mids[which(res$counts == max(res$counts))], sep = ": ")
  } else if (fitgood1[2,1] == "exponential"){
    xsq2 = chi.exp$statistic
    df2 = chi.exp$parameter
    p2 = chi.exp$p.value
    d2 = ks.exp$statistic
    ksp2 = ks.exp$p.value
    xparam21 = paste("Rate", ft.expo$estimate[1], sep = ": ")
  } else if (fitgood1[2,1] == "poisson"){
    xsq2 = chi.pois$statistic
    df2 = chi.pois$parameter
    p2 = chi.pois$p.value
    d2 = ks.pois$statistic
    ksp2 = ks.pois$p.value
    xparam21 = paste("Rate", ft.pois$estimate[1], sep = ": ")
  } else if (fitgood1[2,1] == "lognormal"){
    xsq2 = chi.lognorm$statistic
    df2 = chi.lognorm$parameter
    p2 = chi.lognorm$p.value
    d2 = ks.lognorm$statistic
    ksp2 = ks.lognorm$p.value
    xparam21 = paste("Mean", ft.lognorm$estimate[1], sep = ": ")
    xparam22 = paste("SD", ft.lognorm$estimate[2], sep = ": ")
  } else if (fitgood1[2,1] == "weibull"){
    xsq2 = chi.wei$statistic
    df2 = chi.wei$parameter
    p2 = chi.wei$p.value
    d2 = ks.wei$statistic
    ksp2 = ks.wei$p.value
    xparam21 = paste("Scale", ft.wei$estimate[1], sep = ": ")
    xparam22 = paste("Shape", ft.wei$estimate[2], sep = ": ")
  } else if (fitgood1[2,1] == "gamma"){
    xsq2 = chi.gamma$statistic
    df2 = chi.gamma$parameter
    p2 = chi.gamma$p.value
    d2 = ks.gamma$statistic
    ksp2 = ks.gamma$p.value
    xparam21 = paste("Shape", ft.gamma$estimate[1], sep = ": ")
    xparam22 = paste("Scale", ft.gamma$estimate[2], sep = ": ")
  } else if (fitgood1[2,1] == "uniform"){
    xsq2 = chi.uni$statistic
    df2 = chi.uni$parameter
    p2 = chi.uni$p.value
    d2 = ks.uni$statistic
    ksp2 = ks.uni$p.value
    xparam21 = paste("Min", min(x), sep = ": ")
    xparam22 = paste("Max", max(x), sep = ": ")
  }
  
  if(fitgood1[3,1] == "normal"){
    xsq3 = chi.norm$statistic
    df3 = chi.norm$parameter
    p3 = chi.norm$p.value
    d3 = ks.norm$statistic
    ksp3 = ks.norm$p.value
    xparam31 = paste("Mean", ft.norm$estimate[1], sep = ": ")
    xparam32 = paste("SD", ft.norm$estimate[2], sep = ": ")
  } else if (fitgood1[3,1] == "triangular"){
    xsq3 = chi.triang$statistic
    df3 = chi.triang$parameter
    p3 = chi.triang$p.value
    d3 = ks.triang$statistic
    ksp3 = ks.triang$p.value
    
    xparam31 = paste("Min", min(x), sep = ": ")
    xparam32 = paste("Max", max(x), sep = ": ")
    xparam33 = paste("Mode", res$mids[which(res$counts == max(res$counts))], sep = ": ")
  } else if (fitgood1[3,1] == "exponential"){
    xsq3 = chi.exp$statistic
    df3 = chi.exp$parameter
    p3 = chi.exp$p.value
    d3 = ks.exp$statistic
    ksp3 = ks.exp$p.value
    xparam31 = paste("Rate", ft.expo$estimate[1], sep = ": ")
  } else if (fitgood1[3,1] == "poisson"){
    xsq3 = chi.pois$statistic
    df3 = chi.pois$parameter
    p3 = chi.pois$p.value
    d3 = ks.pois$statistic
    ksp3 = ks.pois$p.value
    xparam31 = paste("Rate", ft.pois$estimate[1], sep = ": ")
  } else if (fitgood1[3,1] == "lognormal"){
    xsq3 = chi.lognorm$statistic
    df3 = chi.lognorm$parameter
    p3 = chi.lognorm$p.value
    d3 = ks.lognorm$statistic
    ksp3 = ks.lognorm$p.value
    xparam31 = paste("Mean", ft.lognorm$estimate[1], sep = ": ")
    xparam32 = paste("SD", ft.lognorm$estimate[2], sep = ": ")
  } else if (fitgood1[3,1] == "weibull"){
    xsq3 = chi.wei$statistic
    df3 = chi.wei$parameter
    p3 = chi.wei$p.value
    d3 = ks.wei$statistic
    ksp3 = ks.wei$p.value
    xparam31 = paste("Scale",ft.wei$estimate[1], sep = ": ")
    xparam32 = paste("Shape", ft.wei$estimate[2], sep = ": ")
  } else if (fitgood1[3,1] == "gamma"){
    xsq3 = chi.gamma$statistic
    df3 = chi.gamma$parameter
    p3 = chi.gamma$p.value
    d3 = ks.gamma$statistic
    ksp3 = ks.gamma$p.value
    xparam31 = paste("Shape", ft.gamma$estimate[1], sep = ": ")
    xparam32 = paste("Scale", ft.gamma$estimate[2], sep = ": ")
  } else if (fitgood1[3,1] == "uniform"){
    xsq3 = chi.uni$statistic
    df3 = chi.uni$parameter
    p3 = chi.uni$p.value
    d3 = ks.uni$statistic
    ksp3 = ks.uni$p.value
    xparam31 = paste("Min", min(x), sep = ": ")
    xparam32 = paste("Max", max(x), sep = ": ")
  }
  
  ########################
  ui <- dashboardPage(
    dashboardHeader(title = "Fitted distributions"),
    dashboardSidebar(disable = TRUE),
    
    dashboardBody(
      box(title = "Table of distributions",width = 6,
          column(width = 6, "Distribution:"),
          column(width = 6, "Square error:"),br(),
          
          column(width = 6, fitgood1[1,1]),
          column(width = 6, fitgood1[1,2]),br(),
          
          column(width = 6, fitgood1[2,1]),
          column(width = 6, fitgood1[2,2]),br(),
          
          column(width = 6, fitgood1[3,1]),
          column(width = 6, fitgood1[3,2]),br(),
          
          column(width = 6, fitgood1[4,1]),
          column(width = 6, fitgood1[4,2]),br(),
          
          column(width = 6, fitgood1[5,1]),
          column(width = 6, fitgood1[5,2]),br(),
          
          column(width = 6, fitgood1[6,1]),
          column(width = 6, fitgood1[6,2]),br(),
          
          column(width = 6, fitgood1[7,1]),
          column(width = 6, fitgood1[7,2]),br(),
          
          column(width = 6, fitgood1[8,1]),
          column(width = 6, fitgood1[8,2])),
      fluidRow(
        tabBox(
          title = "Tests of distributions",
          # The id lets us use input$tabset1 on the server to find the current tab
          id = "tabset1", height = "220px",
          tabPanel("Best fit", plotOutput("plot1",height=300),br(),"Chi-squared test",br(),  "X-squared: ", xsq1,br(),"df: ", df1,br(),"P-value: ",p1,br(),br(),"Kolmogorov-Smirnov Test",br(), "Test statistic: ", d1,br(),"P-value: ", ksp1,br(),br(),"Parameters for the distribution: ", br(), xparam11,br(),xparam12,br(),xparam13),
          tabPanel("2nd Best Fit", plotOutput("plot2",height=300),br(),"Chi-squared test",br(),"X-squared: ", xsq2,br(),"df: ", df2,br(),"P-value: ",p2,br(),br(),"Kolmogorov-Smirnov Test",br(), "Test statistic: ", d2,br(),"P-value: ", ksp2,br(),br(),"Parameters for the distribution: ", br(), xparam21,br(),xparam22,br(),xparam33),
          tabPanel("3rd Best Fit", plotOutput("plot3",height=300),br(),"Chi-squared test",br(),"X-squared: ", xsq3,br(),"df: ", df3,br(),"P-value: ",p3,br(),br(),"Kolmogorov-Smirnov Test",br(), "Test statistic: ", d3,br(),"P-value: ", ksp3,br(),br(),"Parameters for the distribution: ", br(), xparam31,br(),xparam32,br(),xparam33)
        )
      ),
      fluidRow(
        tabBox(
          title = "ECDF",
          id = "tabset2",
          tabPanel("Best fit",plotOutput("plotECDF1",height = 300)),
          tabPanel("2nd Best Fit",plotOutput("plotECDF2",height = 300)),
          tabPanel("3rd Best Fit",plotOutput("plotECDF3",height = 300))
          
        )
      )
    ),
    skin = "green"
  )
  
  server <- function(input,output) {
    
    output$plot1 <- renderPlot({
      hist(dat1, density=20, breaks=seq(xmin,xmax,l=bins+1), prob=TRUE, xlab="x-variable", main="curve over histogram")
      if(fitgood1[1,1] == "normal"){
        curve(dnorm(x,ft.norm$estimate[1], ft.norm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "triangular"){
        curve(dtriangle(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "exponential"){
        curve(dexp(x,ft.expo$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "poisson"){
        curve(dpois(x,ft.pois$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "lognormal"){
        curve(dlnorm(x,ft.lognorm$estimate[1],ft.lognorm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "weibull"){
        curve(dweibull(x,ft.wei$estimate[1],ft.wei$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "gamma"){
        curve(dgamma(x,ft.gamma$estimate[1],ft.gamma$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "uniform"){
        curve(dexp(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      }
    })
    output$plot2 <- renderPlot({
      hist(dat1, density=20, breaks=seq(xmin,xmax,l=bins+1), prob=TRUE, xlab="x-variable", main="curve over histogram")
      if(fitgood1[2,1] == "normal"){
        curve(dnorm(x,ft.norm$estimate[1], ft.norm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "triangular"){
        curve(dtriangle(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "exponential"){
        curve(dexp(x,ft.expo$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "poisson"){
        curve(dpois(x,ft.pois$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "lognormal"){
        curve(dlnorm(x,ft.lognorm$estimate[1],ft.lognorm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "weibull"){
        curve(dweibull(x,ft.wei$estimate[1],ft.wei$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "gamma"){
        curve(dgamma(x,ft.gamma$estimate[1],ft.gamma$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "uniform"){
        curve(dexp(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      }
    })
    output$plot3 <- renderPlot({
      hist(dat1, density=20, breaks=seq(xmin,xmax,l=bins+1), prob=TRUE, xlab="x-variable", main="curve over histogram")
      if(fitgood1[3,1] == "normal"){
        curve(dnorm(x,ft.norm$estimate[1], ft.norm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "triangular"){
        curve(dtriangle(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "exponential"){
        curve(dexp(x,ft.expo$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "poisson"){
        curve(dpois(x,ft.pois$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "lognormal"){
        curve(dlnorm(x,ft.lognorm$estimate[1],ft.lognorm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "weibull"){
        curve(dweibull(x,ft.wei$estimate[1],ft.wei$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "gamma"){
        curve(dgamma(x,ft.gamma$estimate[1],ft.gamma$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "uniform"){
        curve(dexp(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      }
    })
    output$plotECDF1 <- renderPlot({
      plot(ecdf(x), main="ECDF of the data")
      if(fitgood1[1,1] == "normal"){
        curve(pnorm(x,ft.norm$estimate[1], ft.norm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "triangular"){
        curve(ptriangle(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "exponential"){
        curve(pexp(x,ft.expo$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "poisson"){
        curve(ppois(x,ft.pois$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "lognormal"){
        curve(plnorm(x,ft.lognorm$estimate[1],ft.lognorm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "weibull"){
        curve(pweibull(x,ft.wei$estimate[1],ft.wei$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "gamma"){
        curve(pgamma(x,ft.gamma$estimate[1],ft.gamma$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[1,1] == "uniform"){
        curve(pexp(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      }
    })
    output$plotECDF2 <- renderPlot({
      plot(ecdf(x), main="ECDF of the data")
      if(fitgood1[2,1] == "normal"){
        curve(pnorm(x,ft.norm$estimate[1], ft.norm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "triangular"){
        curve(ptriangle(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "exponential"){
        curve(pexp(x,ft.expo$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "poisson"){
        curve(ppois(x,ft.pois$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "lognormal"){
        curve(plnorm(x,ft.lognorm$estimate[1],ft.lognorm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "weibull"){
        curve(pweibull(x,ft.wei$estimate[1],ft.wei$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "gamma"){
        curve(pgamma(x,ft.gamma$estimate[1],ft.gamma$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[2,1] == "uniform"){
        curve(pexp(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      }
    })
    output$plotECDF3 <- renderPlot({
      plot(ecdf(x), main="ECDF of the data")
      if(fitgood1[3,1] == "normal"){
        curve(pnorm(x,ft.norm$estimate[1], ft.norm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "triangular"){
        curve(ptriangle(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "exponential"){
        curve(pexp(x,ft.expo$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "poisson"){
        curve(ppois(x,ft.pois$estimate[1]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "lognormal"){
        curve(plnorm(x,ft.lognorm$estimate[1],ft.lognorm$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "weibull"){
        curve(pweibull(x,ft.wei$estimate[1],ft.wei$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "gamma"){
        curve(pgamma(x,ft.gamma$estimate[1],ft.gamma$estimate[2]), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      } else if (fitgood1[3,1] == "uniform"){
        curve(pexp(x,min(x),max(x)), col="darkblue", lwd=2, add=TRUE, yaxt="n")
      }
    })
    
  }
  shinyApp(ui,server)
}