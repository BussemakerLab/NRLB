setwd("~/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/dbsampling/")

#Testing Wang-Landau covergence for MAFK using different types of Markov Chains.
#We want to minimize the number of DeepBind calls (iterations) needed to converge
#to a valid distribution. The Wang-Landau algorithm divides energy space up into
#bins and attempts to estimate the number of states in each bin. To do this, it 
#requires that the Markov process 'visits' each bin roughly the same number of
#times (even if the DOS varies a lot!). One WL convergence criteria requires that 
#the smallest number of visits for any bin (the minimum of the histogram) be 90% 
#of the mean number of visits for all bins (the mean of the histogram). Every time
#this convergence criteria is met, the value of f (an 'annealing' parameter) is
#halved. The WL algorithm ends when f drops below a critical value (which should
# be small). Here, we try to see if "restarting" the chain with a new random seed
#can speed up the process. For this, we use the MAFK DeepBind model with 150 bins
#and test never restarting the chain, restarting after 500, 1000, 2500, and 5000
#iterations. We can display the convergence of the WL algorithm in a plot where
#the ratio of min(histogram)/mean(histogram) is plotted vs. the number of 
#iterations. Every time convergence is reached (ratio >.9), we can add 1 to the 
#ratio. 

#Needed the function below when the data existed
consolidate = function(input) {
  idx = which(is.nan(input$V2))
  lastIdx = nrow(input)
  
  for (currIdx in rev(idx)) {
    input[currIdx:lastIdx,1] = input[currIdx:lastIdx,1] + input[(currIdx-1),1] + 1
    input[currIdx:lastIdx,2] = input[currIdx:lastIdx,2] + 1
    input[currIdx, 2] = 1
  }
  return(input)
}

#The following was needed when the data existed
# zero = consolidate(read.table("zero.tsv", header=FALSE, stringsAsFactors=FALSE))
# five = consolidate(read.table("five.tsv", header=FALSE, stringsAsFactors=FALSE))
# ten = consolidate(read.table("ten.tsv", header=FALSE, stringsAsFactors=FALSE))
# twentyfive = consolidate(read.table("twentyfive.tsv", header=FALSE, stringsAsFactors=FALSE))
# fifty = consolidate(read.table("fifty.tsv", header=FALSE, stringsAsFactors=FALSE))

load("MCSamplingRuns.Rdata")
combined = rbind(cbind(zero, GROUP="-"), cbind(five, GROUP="500"), cbind(ten, GROUP="1000"), cbind(twentyfive, GROUP="2500"), cbind(fifty, GROUP="5000"))
names(combined) = c("Iterations", "Ratio", "RestartInterval")

ggplot(combined, aes(x=Iterations/1000, y=Ratio, color=RestartInterval, lty=RestartInterval)) +
  theme_bw() + 
  geom_line() + 
  xlab("Iterations (in Thousands)") + 
  theme(legend.position=c(.9,.15), legend.title.align = .5)
ggsave(filename = "RatioPlot.pdf", width=12, height=8)