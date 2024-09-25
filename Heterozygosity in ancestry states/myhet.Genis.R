get_het<-function(x) {
het=x[2,2]/(sum(x[,2]))
return(het)}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  for (filename in args) {
    #print(filename)
    dat <- read.table(file = filename, header = FALSE)
    #print(dat)
    myhet <-  get_het(dat)
    cat(paste(gsub("_psmc_genis_AC","",filename),myhet), sep = "\n")
  }
}

main()

