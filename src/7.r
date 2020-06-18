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
find.mysterious.proteins <- function(cDNA){
  output_file <- paste(TPfolder, "Students/proteins.txt", sep = "")
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

requireNamespace("Biostrings", quietly = TRUE)
library("seqinr")

TPfolder <- "~/ECN/Biologie/13-TP-ProjR/"
mysterious_filename <- paste(TPfolder, "Students/Mysterious_seq.txt", sep = "")
codon_filename <- paste(TPfolder, "Code_shortcut/Genetic_code.txt", sep = "")
reference_filename <- paste(TPfolder, 'Students/PKD1-PKD2_seq.txt', sep = "")
proteins_filename <- paste(TPfolder, 'Students/Protein_sequences.txt', sep = "")
alignment_filename <- paste(TPfolder, 'Students/alignment.txt', sep = "")

ref.prot.sequences <- Biostrings::readAAStringSet(reference_filename)
target.prot.sequences <- Biostrings::readAAStringSet(proteins_filename)
myst_sequences <- read.fasta(file = mysterious_filename, seqtype = "DNA", forceDNAtolower = T)

myst_sequence <- myst_sequences$seq0
myst_proteins <- find.mysterious.proteins(myst_sequence)
biggest_myst_protein <- myst_proteins[which.max(nchar(myst_proteins))]
pair.alignment = Biostrings::pairwiseAlignment(pattern = ref.prot.sequences,
                                               subject = biggest_myst_protein, 
                                               substitutionMatrix = "BLOSUM62",
                                               type = "global")
if( which.max(pair.alignment@score) == 1 ){
  print("[ PKD1 - 4303aa ] is the sequence")
}else{
  print("[ PKD2 - 968aa ] is the sequence")
}


mysequence <- target.prot.sequences$seq12
pair.alignment = Biostrings::pairwiseAlignment(pattern = ref.prot.sequences,
                                               subject = mysequence, 
                                               substitutionMatrix = "BLOSUM62",
                                               type = "local")
Biostrings::writePairwiseAlignments(pair.alignment, alignment_filename)
if( which.max(pair.alignment@score) == 1 ){
  print("[ PKD1 - 4303aa ] is the sequence")
}else{
  print("[ PKD2 - 968aa ] is the sequence")
}




len = length(target.prot.sequences)
matrix.scores = matrix(0, len, len)
for(i in 1:len){
  for(j in 1:len){
    pair.alignment <- Biostrings::pairwiseAlignment(pattern = target.prot.sequences[i],
                                                    subject = target.prot.sequences[j],
                                                    substitutionMatrix = "BLOSUM62",
                                                    type = "global")
    matrix.scores[i, j] = pair.alignment@score
  }
}
heatmap(matrix.scores, symm = T, main = "Alignment of the scores of 32 proteins")


if(FALSE){
output$align.out <- renderPrint({
  ref.sequence <- clean.sequence(ref.prot.sequences[input$select.ref.seq])
  in.seq <- clean.sequence(input$subject.seq)
  if(is.null(in.seq))
    return(NULL)
  pair.alignment <- pairwiseAlignment(attern = ref.sequence,
                                     subject = in.seq,
                                     substitutionMatrix = input$align.matrix,
                                     type = input$align.type)
  return(pair.alignment)
})
}

if(FALSE){
tabPanel("Sequence alignment",
         br(),
         sidebarLayout(sidebarPanel(
           textAreaInput(inputId = "subject.seq",value = "MARVPRVRPPHGFALFLAKEEARKVKRLHGMLRSLLVYMLFLLVTLLASYGDASCHGHAYRLQSAIKQELH",
                         label = "Sequence to align against PKD1 of PKD2 references",
                         resize = 'both',
                         cols = 80,
                         rows = 10),
           selectInput("select.ref.seq", "Select reference sequence", c("PKD1_protein", "PKD2_protein"), "PKD1_protein", multiple = F),
           radioButtons("align.type", "Alignment type", c("global", "local"), "global"),
           selectInput("align.matrix", "Substitution matrix", c("BLOSUM45", "BLOSUM50",  "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250"), "BLOSUM62"),
           sliderInput("slider1", "LabelSlider1", 0, 100, 10),
           sliderInput("slider2", "labelSlider2", 0, 100, 4)
         ),
         mainPanel( h3("Alignment result"), verbatimTextOutput("align.out") )
         )
)
}