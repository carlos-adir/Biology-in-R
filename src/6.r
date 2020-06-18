complementary <- function(seq){
  bases <- "acgt"
  bases_rev <- "tgca"
  complement_seq <- chartr(bases, bases_rev, seq)
  return(rev(complement_seq)) # To adjust from 5->3 again 
}



dna2peptide <- function(sequence, frame=0, sens='F'){
  codon.table <- read.table(file = codon_filename, col.names = c("codon", "aa", "letter"),
                            stringsAsFactors = F)
  codon.table$codon <- tolower(codon.table$codon) # use only lower case characters
  rownames(codon.table) <- codon.table$codon # use codons as table keys

  if(sens == 'R'){
    sequence <- complementary(sequence)
  }
  step = 3
  len <- length(sequence)
  l <- step*((len - frame)%/%step)
  init = frame + 1
  end = l + frame
  seq_coupe <- split(sequence[init:end], ceiling(seq(l)/step))
  for( i in 1:length(seq_coupe) ){
    seq_coupe[i] <- paste(unlist(seq_coupe[i]), collapse='')
  }
  peptides <- codon.table[unlist(seq_coupe), "letter"]
  return( paste(peptides, collapse = "") )
}
sequence <- c('a', 't', 'g', 't', 't', 'c', 't', 't', 't', 'a')
dna2peptide(sequence, 0, 'F')

find.mysterious.proteins <- function(cDNA){
  output_file <- "~/ECN/Biologie/13-TP-ProjR/Students/proteins.txt"
  pattern = 'M[ACDEFGHIKLMNPQRSTVWY]{80,}(X|$)'
  
  proteins <- c()
  step <- 3
  for( frame in 0:(step-1) ){
    for( sens in c("L", "R") ){
      supposed_peps <- paste(dna2peptide(cDNA, frame, sens), collapse = "")
      reg <- gregexpr(pattern, supposed_peps)
      match <- regmatches(supposed_peps, reg)[[1]]
      proteins <- c(proteins, match)
    }
  }
  write(proteins, file = output_file)
  return(proteins)
}

library("seqinr")
codon_filename <- "~/ECN/Biologie/13-TP-ProjR/Code_shortcut/Genetic_code.txt"
filename <- "~/ECN/Biologie/13-TP-ProjR/Students/Mysterious_seq.txt"
sequences <- read.fasta(file = filename, seqtype = "DNA", forceDNAtolower = T)
mysequence <- sequences$seq0
proteins <- find.mysterious.proteins(mysequence)
print(proteins)

output$aa.sequences <- renderUI({
  inFile.2 <- input$dna.file.2
  if (is.null(inFile.2))
    return(NULL)
  dna.sequences <- read.fasta(file = inFile.2$datapath,
                              seqtype = "DNA", forceDNAtolower = T)
  list.of.proteins <- find.mysterious.proteins(dna.sequences$seq0)
  HTML(proteins.to.html(list.of.proteins))
})