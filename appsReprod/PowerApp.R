
library("shiny")
library("shinyWidgets")
library("sampleSize")
source("helpers_power.R")

# Define User Interface ----
ui <- fluidPage(
  titlePanel("Replication Probabilities"),
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(6, 
               numericInput("c", 
                            h4("Relative sample size"), 
                            # value = NA, min=0, max=1)),
                             value = 1, min=0, max=50, step=0.2)),
        
        column(6, 
               numericInput("alpha", 
                            h4("Significance level"), 
                             value = 0.05, min=0, max=1, step=0.01))),
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
      plotOutput("powerPlot", click = "data_click"),
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
  

  output$powerPlot <- renderPlot({ 
    
     # values <- reactiveValues(x=0.05, y=50)
    p.min <- 0.00001
    p.max <- input$p.max
    powerPlot(p.max=p.max)
    p <- exp(seq(log(p.min), log(p.max), length.out=250))
    

    mycol <- c(4,3)
    # if both priors are selected
    if(length(input$prior) >1){
      legend("bottomright", c("Uncertainty No", "Uncertainty Yes"), col=mycol, lty=1)
      # legend("topright", c("Point prior", "Normal prior"), col=mycol, lty=1)
      # standardPowerCurve(p=p, c=input$c, alpha=input$alpha)
      # classicalPowerCurve(p=p, c=input$c, alpha=input$alpha)
    }
     # else{
      if(1 %in% input$prior){
        standardPowerCurve(p=p, c=input$c, alpha=input$alpha, col=mycol[1])
         # reactive(
           segm.power(p=input$p.o, c=input$c, alpha=input$alpha, prior="point")
        points(x=input$p.o, y=100*Power(p=input$p.o, c=input$c, alpha=input$alpha, prior="point"), 
               col=mycol[1], pch=19)
        # )
      }
       if(2 %in% input$prior){
       classicalPowerCurve(p=p, c=input$c, alpha=input$alpha, col=mycol[2])
       segm.power(p=input$p.o, c=input$c, alpha=input$alpha, prior="normal")
       points(x=input$p.o, y=100*Power(p=input$p.o, c=input$c, alpha=input$alpha, prior="normal"), 
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
        power <- 100*StandardPower(n=50, delta=p.to.z(input$p.o)/sqrt(50/(input$c)), sd=1,
                               sig.level=input$alpha, type="one.sample")
        power.r <- round(power, 1)
        return(paste0("Standard Power: ", power.r, "%"))
        }
#         if(2 %in% input$prior){
#           
#         }
#         )
    })
  
  output$info2 <- renderText({
    if(2 %in% input$prior){
              power <- 100*ClassicalPower(n0=50, n=input$c*50, 
                                          delta=p.to.z(input$p.o)/sqrt(50), sd=1,
                                     sig.level=input$alpha, type="one.sample")
              power.r <- round(power, 1)
              return(paste0("Power under uncertainty: ", power.r, "%"))
              
    }
    #         }
    #         )
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
