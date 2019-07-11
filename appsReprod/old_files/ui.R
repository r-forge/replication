
library("shiny")
library("shinyWidgets")
# library("sampleSize")
# library("pCalibrate")


# Define User Interface ----
ui <- fluidPage(
  titlePanel("Design of Replication Studies"),
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(12,
               radioButtons("fixed",
                            h4("Specify"),
                            choices = list("Relative sample size" = 1,
                                           "Replication power" = 2),
                            # "Power" = 2),
                            selected = 1))
      ),
      fluidRow(
        column(6, 
               conditionalPanel(
                 condition = "input.fixed == 1",
                 numericInput("c", 
                              h4("Relative sample size"), 
                              value = 1, min=0, max=50, step=0.2)),
               # column(6, 
               conditionalPanel(
                 condition = "input.fixed == 2",
                 numericInput("power", 
                              h4("Replication power (in %)"),
                              # h4("Power (in %)"), 
                              value = 80, min=0, max=100, step=5))),
        column(6, 
               numericInput("alpha", 
                            h4("Significance level"), 
                            value = 0.05, min=0, max=1, step=0.01))
      ),
      
      fluidRow(
        column(12,
               numericInput("p.o",
                            h4("p-value in original study"),
                            # h4("p-value of original study"),
                            value = 0.01, min=0, max=1, step=0.01))
      ),
      fluidRow(
        column(12, 
               pickerInput("prior", 
                           h4("Power calculation"),
                           # h4("Assumptions about effect size"),
                           # h4("Uncertainty of effect size incorporated"), 
                           choices = list("Traditional" = 1, 
                                          "Predictive" = 2), 
                           selected = 1, multiple=TRUE))),
      
      
      fluidRow(
        column(12,
               sliderInput("p.max",
                           label=h4("Maximum p-value (x-axis)"),
                           value = 0.05, min=0.01, max=0.3, step=0.01))
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
      width=7
    )
  )
)