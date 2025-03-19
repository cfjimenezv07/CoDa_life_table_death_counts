library(shiny)
library(tidyverse)
library(leaflet)

Ui <- fluidPage(
  fluidRow(
    column(4,
           wellPanel(
             div(id = "hidden",
                 selectInput(
                   "countrySelector",
                   "Country",
                   c("USA", "France")
                 ),
             ),
             conditionalPanel(
               "input.countrySelector == 'USA'",
               selectInput("stateSelector",
                           "State",
                           sapply(readRDS("names/names_states.rds"), 
                                  function(x) x)
               )
             ),
             conditionalPanel(
               "input.countrySelector == 'France'",
               selectInput("departmentSelector",
                           "Department",
                           sapply(readRDS("names/names_departments.rds"), 
                                  function(x) x)
               )
             ),
             selectInput("genderSelector",
                         "Gender",
                         c("Male", "Female")
             ),
             selectInput(
               "pcaSelector",
               "Selection method for the FPCA",
               c("EVR", "K = 6")
             ),
             selectInput(
               "forecastSelector",
               "Forecast horizon",
               c("1", "2","3","4","5","6","7","8","9","10")
             )
           ),
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h3(
               "Forecasting density-valued functional panel data",
               style = "font-weight: bold; text-align: center;"
             )),
           div(
             style = "display: flex; justify-content: center; align-items: center; padding-top:5%",
             h4(
               "Authors",
               br(),
               br(),
               "Cristian F. Jiménez-Varón",
               br(),
               "Department of Mathematics",
               br(),
               "University of York",
               br(),
               br(),
               "Ying Sun",
               br(),
               "CEMSE Division",
               br(),
               "King Abdullah University of Science and Technology",
               br(),
               br(),
               "Han Lin Shang",
               br(),
               "Department of Actuarial Studies and Business Analytics",
               br(),
               "Macquarie University")
           )
    ),
    column(4,
           leafletOutput("map")
    ),
    column(4,
           # Upper part with the table
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h4("Point forecast evaluation", style = "font-weight: bold;")
           ),
           div(DT::dataTableOutput("PFETable")),
           # Bottom part with the plot
           div(
             style = "display: flex; justify-content: center; align-items: center; margin-top: 20px;",
             h4("Holdout data - Point forecast", style = "font-weight: bold;")
           ),
           div(plotOutput("myPlot"))
    )
  ),
  tags$style(type = "text/css",
             ".dataTables_length, .dataTables_filter {display:none;}
             h4 {text-align:center;}
             #map {height: calc(100vh - 30px) !important;}
             td {text-align:center !important;}
             #hidden{display: none}")
)
