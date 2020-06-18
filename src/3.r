library("seqinr")
filename <- "~/ECN/Biologie/13-TP-ProjR/Students/Mysterious_seq.txt"
sequences <- read.fasta(file = filename, seqtype = "DNA", forceDNAtolower = T)
mysequence <- sequences$seq0
print(mysequence[1:50])

print(length(mysequence))

bases <- c("a", "c", "g", "t")
composition <- matrix(integer(4), ncol=4)
for(i in 1:4){
  composition[i] <- sum(mysequence == bases[i])
}
colnames(composition) <- bases
print(as.table(composition))

GCcount = composition[2] + composition[3]
GCcontent = GCcount/length(mysequence)
print(paste("Value GC-content:", GCcontent))

complementary <- function(seq){
  bases <- "acgt"
  bases_rev <- "tgca"
  complement_seq <- chartr(bases, bases_rev, seq)
  return(rev(complement_seq)) # To adjust from 5->3 again 
}

mysequence_complement <- complementary(mysequence)
print("Complementary sequence:")
print(mysequence_complement[1:25])