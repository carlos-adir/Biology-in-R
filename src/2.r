compte <- function(sequence, letter){
	return(sum(sequence == letter))
}

sequence <- c("a","a","t","g","a","g","c","t","a","g","c","t","g")
qtt_a <- compte(sequence, "a")
print(qtt_a)

bases = c("a","c","g","t")
composition = rep(0,4)
for(i in 1:4){
	composition[i] = compte(sequence, bases[i])
}
print(composition)