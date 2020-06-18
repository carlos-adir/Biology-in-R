transcribe <- function(cDNA){
  RNA <- chartr("t", "u", cDNA)
  return(RNA)
}