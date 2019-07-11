library("shiny")
library("shinyWidgets")
library("ReplicationSuccess")
source("helpers_power.R")
source("helpers_relSampleSize.R")

# Define User Interface ----
ui <- fluidPage(
  titlePanel("Planung von Replikationsstudien"), #Design of Replication Studies
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
      column(12,
             radioButtons("fixed",
                                h4("Fixiere"), #Specify
                                choices = list("Relative Stichprobengrösse" = 1, #Relative sample size
                                               "Replikations Power" = 2), #Replication power
                                selected = 1))
      ),
      fluidRow(
        column(6, 
               conditionalPanel(
                 condition = "input.fixed == 1",
                 numericInput("c", 
                              h4("Relative Stichprobengrösse"), #Relative sample size
                              value = 1, min = 0, max = 50, step = 0.2)),
        # column(6, 
               conditionalPanel(
               condition = "input.fixed == 2",
               numericInput("power", 
                            h4("Replikations Power (in %)"), #Replication power (in %)
                            value = 80, min = 0, max = 100, step = 5))),
        column(6, 
               numericInput("alpha", 
                            h4("Signifikanz Niveau"), #Significance level
                            value = 0.05, min = 0, max = 1, step = 0.01))
      ),
        
      fluidRow(
        column(12,
               numericInput("p.o",
                            h4(HTML("<em>p</em>-Wert der Originalstudie (zweiseitig)")), 
                            #p-value in original study
                            value = 0.01, min = 0, max = 1, step = 0.01))
      ),
        fluidRow(
          column(12, 
                 pickerInput("prior", 
                             h4("Power Berechnung"), #Power calculation
                             choices = list("Traditionell" = 1, #Traditional
                                            "Prädiktiv" = 2), #Predictive
                             selected = 1, multiple = TRUE))),
      
    
      fluidRow(
        column(12,
               sliderInput("p.max",
                           label=h4(HTML("Maximaler <em>p</em>-Wert (x-Achse)")), #Maximum p-value (x-axis)
                           value = 0.05, min = 0.01, max = 0.3, step = 0.01))
      )),
    
    mainPanel(
      # You can set many attributes, not just colors
      tags$style(type='text/css', '#info {background-color: white; color: blue;}'),
      tags$style(type='text/css', '#info2 {background-color: white; color: limegreen;}'),
      
      plotOutput("Plot", click = "data_click"),
      # htmlOutput("info"),
      # htmlOutput("info2"),
      verbatimTextOutput("info"),
      verbatimTextOutput("info2"),
      # verbatimTextOutput("warning"),
      htmlOutput("warning"),
      width = 7
    )
  )
)

######################################################################

# Define server logic ----
server <- function(input, output, session) {
  

  ## Plot output
  output$Plot <- renderPlot({ 
    
    p.min <- 0.000001
    p.max <- input$p.max
    
    if(input$fixed == 1){ # if relative sample size fixed
      powerPlot(p.max = p.max)
      p <- exp(seq(log(p.min), log(p.max), length.out = 250))
      mycol <- c(4, 3)
      if(1 %in% input$prior){ # traditional = conditional = point prior
        PowerCurve(p = p, c = input$c, alpha = input$alpha, prior = "conditional",
                   col = mycol[1], lwd = 2)
        segm.power(p = input$p.o, c = input$c, alpha = input$alpha, prior = "conditional")
        points(x = input$p.o, 
               y = 100*powerSignificance(po = input$p.o, c = input$c, level = input$alpha, 
                                         designPrior = "conditional"), col = mycol[1], pch = 19)
      }
      if(2 %in% input$prior){ # predictive = normal prior
        PowerCurve(p = p, c = input$c, alpha = input$alpha, prior = "predictive",
                   col = mycol[2], lwd = 2)
        segm.power(p = input$p.o, c = input$c, alpha = input$alpha, prior = "predictive")
        points(x = input$p.o, 
               y = 100*powerSignificance(po = input$p.o, c = input$c, level = input$alpha,
                                         designPrior = "predictive"), col = mycol[2], pch = 19)
      }
      if(length(input$prior) >1){ # both priors are selected
        legend("topright", c("Traditionell", "Prädiktiv"), col = mycol, 
               lty = 1, lwd = 2, cex = 1.1)
      }
    }
    
    if(input$fixed == 2){ # if power fixed
      if(input$power/100 <= input$alpha)
        stop("Die Power muss grösser sein als das Signifikanz Niveau.")
      #The power has to be larger than the significance level.
      c.max <- 10
      sampleSizePlot(p.max = p.max, c.max = c.max)
      p <- exp(seq(log(p.min), log(p.max), length.out=250))
      mycol <- c(4, 3)
      
        if(1 %in% input$prior){ # traditional = conditional = point prior
          SamSizeCurve(p = p, alpha = input$alpha, power = input$power/100, 
                       prior = "conditional", col = mycol[1], lwd = 2)
          segm.rss(p = input$p.o, power = input$power/100, alpha = input$alpha,
                   prior = "conditional", col = 2)
          points(x = input$p.o, 
                 y = sampleSizeSignificance(po = input$p.o, power = input$power/100,
                                            level = input$alpha, designPrior = "conditional"),
                 col = mycol[1], pch = 19)
        }
         if(2 %in% input$prior){ # predictive = normal prior
           SamSizeCurve(p = p, alpha = input$alpha, power = input$power/100, 
                        prior = "predictive", col = mycol[2], lwd = 2)
           segm.rss(p = input$p.o, power = input$power/100, alpha = input$alpha, 
                    prior = "predictive", col = 2)
           points(x = input$p.o, 
                  y = sampleSizeSignificance(po = input$p.o, power = input$power/100,
                                             level = input$alpha, designPrior = "predictive"), 
                  col = mycol[2], pch = 19)
         }
      if(length(input$prior) >1){ # both priors are selected
        legend("topleft", c("Traditionell", "Prädiktiv"), col = mycol,
               lty = 1, lwd = 2, cex = 1.1)
      }
    }
    
    # update when user clicks on plot
    observeEvent(input$data_click, {
          updateNumericInput(session = session, inputId = "p.o", 
                             value = format(input$data_click$x, digits = 2, nsmall = 1))
        })
  }, height = 400)

  ## Text output under point prior
  output$info <- renderText({
    
    if(1 %in% input$prior){
      if(input$fixed == 1){ # if relative sample size fixed
        power <- 100*powerSignificance(po = input$p.o, c = input$c, level = input$alpha, 
                                       designPrior = "conditional")
        power.r <- round(power, 1)
        return(paste0("Traditionelle Power: ", power.r, "%"))
      }
      if(input$fixed == 2){ # if power fixed
        if(input$power/100 <= input$alpha)
          stop(" ")
        c <- sampleSizeSignificance(po = input$p.o, power = input$power/100, 
                                    level = input$alpha, designPrior = "conditional")
        c.r <- format(c, digits = 2)
        return(paste0("Relative Stichprobengrösse: ", c.r))
      }
    }
  })
  
  ## Text output under normal prior
  output$info2 <- renderText({
    if(2 %in% input$prior){ 
      if(input$fixed == 1){ # if relative sample size fixed
        power <- 100*powerSignificance(po = input$p.o, c = input$c, level = input$alpha, 
                                       designPrior = "predictive")
        power.r <- round(power, 1)
        return(paste0("Prädiktive Power: ", power.r, "%"))
      }
      if(input$fixed == 2){ # if power fixed
        if(input$power/100 <= input$alpha)
          stop(" ")
        c <- sampleSizeSignificance(p = input$p.o, power = input$power/100, 
                                    level = input$alpha, designPrior = "predictive")
        c.r <- format(c, digits = 2)
        return(paste0("Relative Stichprobengrösse: ", c.r))
      }
    }
  })
  
  ## Text output warnings
  output$warning <- renderUI({
    power.max <- 100*(1 - input$p.o/2)
    if(input$fixed == 2 && input$power >= power.max && 2 %in% input$prior){
      power.max.r <- format(power.max, digits = 1, nsmall = 1)
      str1 <- paste0("Warnung: Eine so grosse Power kann nicht erreicht werden.")
      #Warning: Such a high predictive power cannot be attained.
      str2 <- paste0("Die Replikations Power muss kleiner sein als ",
                     power.max.r, "%.")
      #The replication power must be smaller than 
      HTML(paste(str1, str2))
     }
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
