source("mapRender.R")
source("renderPFEtable.R")
# source("renderIFEtable.R")

get.state <- function(input) {
  if(input$countrySelector == "USA") {
    return(input$stateSelector)
  }
  return(input$departmentSelector)
}

get.pca <- function(input) {
  if(input$pcaSelector == "EVR") {
    return(input$pcaSelector)
  } 
  return("K")
}




SERVER <- function(input, output) {
  output$map <- renderLeaflet({
    render.map(input$countrySelector, get.state(input))
  })
  
  output$PFETable <- DT::renderDataTable({
    # Define the custom table header (sketch.upper)
    sketch.upper <- htmltools::withTags(
      table(
        class = 'display',
        thead(
          tr(
            th(colspan = 2, 'FM-ANOVA', style = 'text-align:center;'),
            th(colspan = 2, 'FMP-ANOVA', style = 'text-align:center;'),
            th(colspan = 2, 'TNH', style = 'text-align:center;'),
            th(colspan = 2, 'GSY', style = 'text-align:center;'),
            th(colspan = 2, 'MEM', style = 'text-align:center;')
          ),
          tr(
            lapply(c(rep(c('KLD', 'JSD'), 5)), function(x) th(x, style = 'text-align:center;'))
          )
        )
      )
    )
    
    # Render the table with the custom header
    DT::datatable(
      gen_pfe_table(input$countrySelector, get.state(input), 
                    input$genderSelector, get.pca(input), input$forecastSelector), 
      rownames = FALSE, 
      container = sketch.upper
    )
  })
  
  output$myPlot <- renderPlot({
    data <- gen_pfe_curves_table(input$countrySelector, get.state(input), 
                           input$genderSelector, get.pca(input), input$forecastSelector)
    
    # Plot the first curve in red with lwd=1.5
    plot(1:111, data[,1], type="l", col="red", lwd=2, xlab="Age", ylab="Forecast",ylim=c(min(data),max(data)),lty=1)
  
  # Add subsequent curves in different colors and lwd=1.5
  lines(1:111, data[,2], col="blue", lwd=1.5,lty=2)
  lines(1:111, data[,3], col="purple", lwd=1.5,lty=3)
  lines(1:111, data[,4], col="darkgreen", lwd=1,lty=4)
  lines(1:111, data[,5], col="violet", lwd=1,lty=5)
  lines(1:111, data[,6], col="orange", lwd=1,lty=6)
  
  # Add the legend in the top-left corner
  legend("topleft", legend=c("Holdout data", "FM", "FMP", "TNH", "GSY", "MEM"), 
         col=c("red", "blue", "purple", "darkgreen", "violet", "orange"),lty=c(1, 2, 3, 4, 5, 6), 
         , lwd=1.5,ncol = 2,cex=0.5)  # bty="n" removes the box around the legend
  })

}
