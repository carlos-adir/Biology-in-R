#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(seqinr)
library(Biostrings)

# 'util' file for genetic code logic and strings cleaning
source("util.R")

# server logic to read fasta files
shinyServer(function(input, output) {
  TPfolder <- "~/ECN/Biologie/13-TP-ProjR/"
  reference_filename <- paste(TPfolder, 'Students/PKD1-PKD2_seq.txt', sep = "")
  ref.prot.sequences <- Biostrings::readAAStringSet(reference_filename)
  
  output$dna.content <- DT::renderDataTable({
    
    # save inFile for internal logic
    inFile <- input$dna.file
    
    if (is.null(inFile))
      return(NULL)
    
    dna.sequences <- read.fasta(inFile$datapath, seqtype = "DNA", forceDNAtolower = T)
    
    dna.content <- compute.dna.metrics(dna.sequences)
    
    return(dna.content)
  })
  
  output$aa.sequences <- renderUI({
    inFile.2 <- input$dna.file.2
    if (is.null(inFile.2))
      return(NULL)
    dna.sequences <- read.fasta(file = inFile.2$datapath,
                                seqtype = "DNA", forceDNAtolower = T)
    list.of.proteins <- find.mysterious.proteins(dna.sequences$seq0)
    HTML(proteins.to.html(list.of.proteins))
  })
  
  output$align.out <- renderPrint({
    ref.sequence <- clean.sequence(ref.prot.sequences[input$select.ref.seq])
      in.seq <- clean.sequence(input$subject.seq)
    if(is.null(in.seq))
      return(NULL)
    pair.alignment <- pairwiseAlignment(pattern = ref.sequence,
                                        subject = in.seq,
                                        substitutionMatrix = input$align.matrix,
                                        type = input$align.type)
    
      return(pair.alignment)
  })
  
})
