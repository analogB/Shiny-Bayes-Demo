library(shiny)
require(rjags)
source("DBDA2E-utilities.R")

ui <- fluidPage( headerPanel('Bayesian Analysis of Sample Proportion'),
                 sidebarPanel(
                 numericInput(n, 'Sample Size', 100, min = 1, max = NA, step = 1,width = NULL),
                 numericInput(nTrue, 'Observed Positives', 50, min = 1, max = n, step = 1, width = NULL),
                 renderText({c(n-nTrue, ' Observed Negatives')}),
                 numericInput(CI, 'Confidence Interval (2-sided)', 0.9, min = 0.5, max = 0.99, step = 0.01, width = NULL)
                 ),
                 mainPanel(
                   plotOutput('diagnostics'),
                   renderText({c(lowerEstimate, ' Lower Estimate')}),
                   renderText({c(meanEstimate,  ' Mean Estimate' )}),
                   renderText({c(upperEstimate, ' Upper Estimate')})
                 ))

server <- function(input, output) {}

shinyApp(ui = ui, server = server)


CI <- 0.9
n <- 50
nTrue <- 35
nFalse <- n-nTrue


# FORMAT DATA
myData = c(rep(0,nFalse),rep(1,nTrue))
y <- myData
nTotal <- length(y)

#Data list for rjags
dataList = list(
  y = y,
  nTotal = nTotal
)

#------------------------------------------------------------------------------
# RUN CHAINS
# initsList = list( theta=c(0.5,0.5,0.5) , m=3 ) #not necessary with simple models
parameters    = c("theta") 
adaptSteps    =  1000      # Number of steps to "tune" the samplers.
burnInSteps   = 10000      # Number of steps to "burn-in" the samplers.
numSavedSteps =  2000      # Total number of steps to save collectively from all chains. 
nChains       =     3      # Number of chains to run.
thinSteps     =     1      # Number of steps to "thin" (1=keep every step).
nPerChain     = ceiling((numSavedSteps*thinSteps)/nChains) # Steps per chain.

# Create, initialize, and adapt the model:
jagsModel = jags.model( "proportion_model.txt" , data=dataList ,  n.chains=nChains , n.adapt=adaptSteps )  #inits=initsList ,
cat("Burning in the MCMC chain...\n") # Burn-in:
update(jagsModel,n.iter=burnInSteps)
cat("Sampling final MCMC chain...\n") # The saved MCMC chain:
codaSamples <- coda.samples(jagsModel,variable.names=parameters,n.iter=nPerChain,thin=thinSteps)
diagnostics <- diagMCMC(codaSamples, parName=paste('theta'))

mcmcMat = as.matrix(codaSamples,chains=TRUE)
summaryInfo = summarizePost( mcmcMat[,"theta"],credMass=CI)
meanEstimate <- summaryInfo[grep('Mean',names(summaryInfo))] 
upperEstimate <- summaryInfo[grep('HDIhigh',names(summaryInfo))] 
lowerEstimate <- summaryInfo[grep('HDIlow',names(summaryInfo))] 



