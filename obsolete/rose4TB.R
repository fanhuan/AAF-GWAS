# set the root sequence
filenamesequence <- "/Users/iveslab/Documents/projects/GWAS/TB/Mcanettii.fasta"
rootgenome <- scan(filenamesequence, what="character", sep="", quiet=TRUE)

# This starts to construct the Rosefile from a universal header
system("cp /Users/iveslab/Documents/projects/GWAS/TB/Rosefile_TB_header /Users/iveslab/Dropbox/RESEARCH/AAF/Code_Sources/rose-1.3/Rosefile_TB_20141125")

# This adds the root sequence into Rosefile_TB
filenametemp <- "/Users/iveslab/Dropbox/RESEARCH/AAF/Code_Sources/rose-1.3/Rosefile_TB_20141125"
cat(rootgenome, file=filenametemp, append=T)
cat("\"", file=filenametemp, append=T, sep="")

#rose simulation
rose.call <- "cd /Users/iveslab/Dropbox/RESEARCH/AAF/Code_Sources/rose-1.3; ./rose Rosefile_TB_20141125"
system(rose.call, ignore.stdout = T)
simulated_seqs <- scan(file = "/Users/iveslab/Dropbox/RESEARCH/AAF/Code_Sources/rose-1.3/simfiles.fas", what = "character", quiet=T)
dataindex <- grep(">",simulated_seqs)

# construct file for phylokmer
nspp <- length(dataindex)
dataindex <- c(dataindex,length(rosefile))
system("rm /Users/iveslab/Documents/projects/GWAS/TB/RoseSim/*")
for(i in 1:nspp){
  spname <- substr(rosefile[dataindex[i]],start=2,stop=nchar(rosefile[dataindex[i]]))
  sequence <- paste(rosefile[(dataindex[i]+1):(dataindex[i+1]-1)], collapse="")
  
  filename1 <- paste(c("/Users/iveslab/Documents/projects/GWAS/TB/RoseSim/", spname, ".fa"), collapse="")
  write(rosefile[dataindex[i]], file=filename1)
  write(sequence, file=filename1, append=T)
}
