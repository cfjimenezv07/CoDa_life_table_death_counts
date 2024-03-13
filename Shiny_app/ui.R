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
               "coverageSelector",
               "Nominal Coverage",
               c("80%", "95%")
             ),
             selectInput(
               "pcaSelector",
               "Selection method for the FPCA",
               c("EVR", "K = 6")
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
               "Cristian F. Jiménez-Varón and Ying Sun",
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
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h4("Point forecast evaluation", style = "font-weight: bold;")
           ),
           div(DT::dataTableOutput("PFETable")),
           # DT::dataTableOutput("summaryMetricsTable"),
           conditionalPanel(
             "input.coverageSelector == '80%'",
             div(
               style = "display: flex; justify-content: center; align-items: center",
               h4("Interval forecast evaluation 80% nominal coverage",
                  style = "font-weight: bold;")
             )
           ),
           conditionalPanel(
             "input.coverageSelector == '95%'",
             div(
               style = "display: flex; justify-content: center; align-items: center",
               h4("Interval forecast evaluation 95% nominal coverage",
                  style = "font-weight: bold;")
             )
           ),
           div(DT::dataTableOutput("IFETable")),
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h6("* All values were multiplied by a factor of 100.",
                style = "font-weight: bold;")
           )
    )
  ),
  tags$style(type = "text/css",
             ".dataTables_length, .dataTables_filter {display:none;}
             h4 {text-align:center;}
             #map {height: calc(100vh - 30px) !important;}
             td {text-align:center !important;}
             #hidden{display: none")
)