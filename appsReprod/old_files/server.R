
source("helpers_power.R")
source("helpers_relSampleSize.R")


# Define server logic ----
server <- function(input, output, session) {
  
  
  output$Plot <- renderPlot({ 
    
    p.min <- 0.000001
    p.max <- input$p.max
    
    if(input$fixed == 1){
      
      powerPlot(p.max=p.max)
      p <- exp(seq(log(p.min), log(p.max), length.out=250))
      
      
      mycol <- c(4,3)
      
      if(1 %in% input$prior){
        standardPowerCurve(p=p, c=input$c, alpha=input$alpha, col=mycol[1])
        segm.power(p=input$p.o, c=input$c, alpha=input$alpha, prior="point")
        points(x=input$p.o, y=100*Power(p=input$p.o, c=input$c, alpha=input$alpha, prior="point"), 
               col=mycol[1], pch=19)
      }
      if(2 %in% input$prior){
        classicalPowerCurve(p=p, c=input$c, alpha=input$alpha, col=mycol[2])
        segm.power(p=input$p.o, c=input$c, alpha=input$alpha, prior="normal")
        points(x=input$p.o, y=100*Power(p=input$p.o, c=input$c, alpha=input$alpha, prior="normal"), 
               col=mycol[2], pch=19)
      }
      # if both priors are selected
      if(length(input$prior) >1){
        legend("bottomright", c("Traditional", "Predictive"), col=mycol, lty=1, lwd=2, cex=1.1)
      }
    }
    if(input$fixed == 2){
      if(input$power/100 <= input$alpha)
        stop("The power has to be larger than the significance level.")
      c.max <- 10
      sampleSizePlot(p.max=p.max, c.max=c.max)
      p <- exp(seq(log(p.min), log(p.max), length.out=250))
      
      
      mycol <- c(4,3)
      if(1 %in% input$prior){
        SamSizeCurve(p=p, alpha=input$alpha, power=input$power/100, prior="point", col=mycol[1])
        segm.rss(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="point", col=2)
        points(x=input$p.o, y=relSampleSize(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="point"), 
               col=mycol[1], pch=19)
      }
      if(2 %in% input$prior){
        SamSizeCurve(p=p, alpha=input$alpha, power=input$power/100, prior="normal", col=mycol[2])
        segm.rss(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="normal", col=2)
        points(x=input$p.o, y=relSampleSize(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="normal"), 
               col=mycol[2], pch=19)
      }
      # if both priors are selected
      if(length(input$prior) >1){
        legend("topleft", c("Traditional", "Predictive"), col=mycol, lty=1, lwd=2, cex = 1.1)
      }
    }
    
    
    
    observeEvent(input$data_click, {
      updateNumericInput(session, "p.o", value= format(input$data_click$x, digits=2, nsmall=1))
    })
  }, height=400)
  #   
  output$info <- renderText({
    
    if(1 %in% input$prior){
      if(input$fixed == 1){
        power <- 100*StandardPower(n=50, delta=p.to.z(input$p.o)/sqrt(50/(input$c)), sd=1,
                                   sig.level=input$alpha, type="one.sample")
        power.r <- round(power, 1)
        return(paste0("Traditional power: ", power.r, "%"))
      }
      if(input$fixed == 2){
        if(input$power/100 <= input$alpha)
          stop(" ")
        c <- relSampleSize(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="point")
        c.r <- format(c, digits=2)
        return(paste0("Relative sample size: ", c.r))
        # return(paste0("Relative sample size ignoring uncertainty: ", c.r))
      }
    }
  })
  
  output$info2 <- renderText({
    if(2 %in% input$prior){
      if(input$fixed == 1){
        power <- 100*ClassicalPower(n0=50, n=input$c*50, 
                                    delta=p.to.z(input$p.o)/sqrt(50), sd=1,
                                    sig.level=input$alpha, type="one.sample")
        power.r <- round(power, 1)
        return(paste0("Predictive power: ", power.r, "%"))
        # return(paste0("Power under uncertainty: ", power.r, "%"))
      }
      if(input$fixed == 2){
        if(input$power/100 <= input$alpha)
          stop(" ")
        c <- relSampleSize(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="normal")
        c.r <- format(c, digits=2)
        return(paste0("Relative sample size: ", c.r))
        # return(paste0("Relative sample size under uncertainty: ", c.r))
      }
    }
  })
  
  output$warning <- renderUI({
    power.max <- 100*(1 - input$p.o/2)
    if(input$fixed==2 && input$power >= power.max && 2 %in% input$prior){
      power.max.r <- format(power.max, digits=1, nsmall=1)
      str1 <- paste0("Warning: Such a high predictive power cannot be attained.")
      str2 <- paste0("The replication power must be smaller than ",
                     power.max.r, "%.")
      HTML(paste(str1, str2))
      # HTML(paste(str1, str2, sep = '<br/>'))
    }
  })
}