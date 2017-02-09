library(shiny)
require(rjags)
source("DBDA2E-utilities.R")

ui <- fluidPage( headerPanel('Bayesian Analysis of Sample Proportion'),
                 sidebarPanel(
                 numericInput('n', 'Sample Size', 100, min = 1, max = NA, step = 1,width = NULL),
                 numericInput('nTrue', 'Observed Positives', 50, min = 1, max = n, step = 1, width = NULL),
                 renderText({c('nFalse', ' Observed Negatives')}),
                 numericInput('CI', 'Confidence Interval (2-sided)', 0.9, min = 0.5, max = 0.99, step = 0.01, width = NULL)
                 ),
                 mainPanel(
                   plotOutput('diagnostics'),
                   renderText({c('lowerEstimate', ' Lower Estimate')}),
                   renderText({c('meanEstimate',  ' Mean Estimate' )}),
                   renderText({c('upperEstimate', ' Upper Estimate')})
                 ))

server <- function(input, output) {
  output$nFalse <- input$n - input$nTrue
  output$y <- c(rep(0,output$nFalse),rep(1,input$nTrue))
  output$nTotal <- length(output$y)
  output$dataList = list(
    y = output$y,
    nTotal = output$nTotal
  )
  parameters    = c("theta") 
  adaptSteps    =  1000      # Number of steps to "tune" the samplers.
  burnInSteps   = 10000      # Number of steps to "burn-in" the samplers.
  numSavedSteps =  2000      # Total number of steps to save collectively from all chains. 
  nChains       =     3      # Number of chains to run.
  thinSteps     =     1      # Number of steps to "thin" (1=keep every step).
  nPerChain     = ceiling((numSavedSteps*thinSteps)/nChains) # Steps per chain.
  
  output$jagsModel = jags.model( "proportion_model.txt" , data=output$dataList ,  n.chains=nChains , n.adapt=adaptSteps )  #inits=initsList ,
  #cat("Burning in the MCMC chain...\n") # Burn-in:
  output$update(output$jagsModel,n.iter=burnInSteps)
  #cat("Sampling final MCMC chain...\n") # The saved MCMC chain:
  output$codaSamples <- coda.samples(output$jagsModel,variable.names=parameters,n.iter=nPerChain,thin=thinSteps)
  output$diagnostics <- diagMCMC(output$codaSamples, parName=paste('theta'))
  
  output$mcmcMat = as.matrix(output$codaSamples,chains=TRUE)
  output$summaryInfo = summarizePost( output$mcmcMat[,"theta"],credMass=input$CI)
  output$meanEstimate <- output$summaryInfo[grep('Mean',names(summaryInfo))] 
  output$upperEstimate <- output$summaryInfo[grep('HDIhigh',names(summaryInfo))] 
  output$lowerEstimate <- output$summaryInfo[grep('HDIlow',names(summaryInfo))] 

}

shinyApp(ui = ui, server = server)