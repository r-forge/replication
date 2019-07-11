
library("shiny")
library("shinyWidgets")
library("sampleSize")
library("pCalibrate")
source("helpers_RelSampleSize.R")

# Define User Interface ----
ui <- fluidPage(
  titlePanel("Relative sample size of the replication study"),
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(6, 
               numericInput("power", 
                            h4("Power (in %)"), 
                            # value = NA, min=0, max=1)),
                            value = 80, min=0, max=100, step=5)),
        
        column(6, 
               numericInput("alpha", 
                            h4("Significance level"), 
                            value = 0.05, min=0, max=1, step=0.01))
      ),
      fluidRow(
        column(12,
               numericInput("p.o",
                            h4("p-value of original study"),
                            # value = NA, min=0, max=1)),
                            value = 0.02, min=0, max=1, step=0.01))
      ),
        fluidRow(
          column(12, 
                 pickerInput("prior", 
                             h4("Uncertainty on effect size incorporated"), 
                             choices = list("No" = 1, 
                                            "Yes" = 2), 
                             selected = 1, multiple=TRUE))),
      # fluidRow(
      #   column(12,
      #          numericInput("p.max",
      #                       h4("Max. p-value"),
      #                       value = 0.05, min=0, max=1, step=0.05))
      # ),
      
      # fluidRow(
      # column(6, 
      #        checkboxGroupInput("prior", 
      #                           h3("Design prior"), 
      #                           choices = list("Point" = 1, 
      #                                          "Normal" = 2), 
      #                           selected = 1)))
      # ),
      fluidRow(
        column(12,
               sliderInput("p.max",
                           label=h4("Max. p-value"),
                           value = 0.05, min=0.01, max=0.3, step=0.01))
      )),
    
    mainPanel(
      plotOutput("SampleSizePlot", click = "data_click"),
      # dblclick = clickOpts("data_dblclick")),
      verbatimTextOutput("info"),
      verbatimTextOutput("info2"),
      width=7
    )
  )
)

######################################################################

# Define server logic ----
server <- function(input, output, session) {
  

  output$SampleSizePlot <- renderPlot({ 
    
     # values <- reactiveValues(x=0.05, y=50)
    p.min <- 0.00001
    p.max <- input$p.max
    c.max <- ceiling(max(relSampleSize(p=p.max, power=input$power/100, alpha=input$alpha, prior="point"),
                 relSampleSize(p=p.max, power=input$power/100, alpha=input$alpha, prior="normal")))
    sampleSizePlot(p.max=p.max, c.max=c.max)
    p <- exp(seq(log(p.min), log(p.max), length.out=250))
    

    mycol <- c(4,3)
    # if both priors are selected
    if(length(input$prior) >1){
      legend("bottomright", c("Uncertainty No", "Uncertainty Yes"), col=mycol, lty=1)
      # legend("topright", c("Point prior", "Normal prior"), col=mycol, lty=1)
    }
     # else{
      if(1 %in% input$prior){
        SamSizeCurve(p=p, alpha=input$alpha, power=input$power/100, prior="point", col=mycol[1])
         # reactive(
           segm.rss(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="point", col=2)
        points(x=input$p.o, y=relSampleSize(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="point"), 
               col=mycol[1], pch=19)
        # )
      }
       if(2 %in% input$prior){
         SamSizeCurve(p=p, alpha=input$alpha, power=input$power/100, prior="normal", col=mycol[2])
         # reactive(
         segm.rss(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="normal", col=2)
         points(x=input$p.o, y=relSampleSize(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="normal"), 
                col=mycol[2], pch=19)
     }
    
        observeEvent(input$data_click, {
          updateNumericInput(session, "p.o", value= format(input$data_click$x, digits=2, nsmall=1))
        })
       #          values$x <- input$data_click$x
       #          # coords$y <- Power(p=coords$x, c=input$c, alpha=input$alpha, prior=input$prior)
       #   # y <- input$data_click$y
       # }
       # )
    #     
  }, height=400)
#   
  output$info <- renderText({
    # mycol <- c(4,3)
    if(1 %in% input$prior){
      c <- relSampleSize(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="point")
      c.r <- format(c, digits=2)
      return(paste0("Relative sample size without uncertainty: ", c.r))
    }
    #         if(2 %in% input$prior){
    #           
    #         }
    #         )
  })
  
  output$info2 <- renderText({
    if(2 %in% input$prior){
      c <- relSampleSize(p=input$p.o, power=input$power/100, alpha=input$alpha, prior="normal")
      c.r <- format(c, digits=2)
      return(paste0("Relative sample size under uncertainty: ", c.r))
      
    }
    #         }
    #         )
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
