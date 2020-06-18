compute.dna.metrics <- function(dna.sequences){
  # Bonus exercise 7.5: combine the 6 sapply calls
  dna.content <- data.frame(names=names(dna.sequences))
  
  dna.content$a <- sapply(dna.sequences, compte, 'a')
  dna.content$c <- sapply(dna.sequences, compte, 'c')
  dna.content$g <- sapply(dna.sequences, compte, 'g')
  dna.content$t <- sapply(dna.sequences, compte, 't')
  
  dna.content$length <- sapply(dna.sequences, length)
  dna.content$GC.content <- sapply(dna.sequences, GC)
  
  return(dna.content)
}

compte <- function(seq, letter){
  count=sum(seq==letter)
  return(count)
}

codon.table <- read.table(file = "../Code_shortcut/Genetic_code.txt", 
                          col.names = c("codon", "aa", "letter"),
                          stringsAsFactors = F)
codon.table$codon <- tolower(codon.table$codon) # use only lower case characters
rownames(codon.table) <- codon.table$codon # use codons as table keys


complementary <- function(seq){
  bases <- "acgt"
  bases_rev <- "tgca"
  complement_seq <- chartr(bases, bases_rev, seq)
  return(rev(complement_seq)) # To adjust from 5->3 again 
}

dna2peptide <- function(sequence, frame=0, sens='F'){
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

clean.sequence <- function(dirty.seq){
  return(gsub('[^A-Z]', replacement = "", toupper(dirty.seq)))
}

proteins.to.html <- function(list.of.proteins){
  # Every 60 amino acids, insert a new line (html code)
  proteins.cut <- gsub("(.{61,}?)", "\\1</br>", list.of.proteins)
  # Use console typo in html
  proteins.in.span <- paste0("<p style='font-family:monospace'>",proteins.cut,"</span>")
  return(paste0(proteins.in.span, collapse = "</br></br>"))
}

# more to come