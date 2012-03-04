## motif.R -- a toy motif finder using Gibbs sampling

library(Biostrings)
NUCLEOTIDES <- c("A", "T", "C", "G")

initPFM <- function(width=8, alphabet=NUCLEOTIDES) {
  matrix(0, nrow=length(alphabet), ncol=width, dimnames=list(alphabet, seq_len(width)))
}

initBackground <- function(sequences, alphabet=NUCLEOTIDES) {
  x <- numeric(length(alphabet))
  names(x) <- alphabet
  x
}

updateBackground <- function(sequences, background) {
  background + colSums(alphabetFrequency(sequences))[names(background)]
}

updatePFM <- function(motif, pfm) {
  stopifnot(ncol(pfm) == nchar(motif))
  letters <- strsplit(motif, '')[[1]]
  for (position in seq_along(letters)) {
    letter <- letters[position]
    pfm[letter, position] <- pfm[letter, position] + 1
  }
  pfm
}

removeFromPWM <- function(motif, pfm) {
  stopifnot(ncol(pfm) == nchar(motif))
  letters <- strsplit(motif, '')[[1]]
  for (position in seq_along(letters)) {
    letter <- letters[position]
    current.count <- pfm[letter, position]
    if (current.count == 0)
      next()
    pfm[letter, position] <- pfm[letter, position] - 1
  }
  pfm
}


samplePositions <- function(sequences, width=8, prob=NULL) {
  starts <- sapply(seq_along(sequences), function(i) {
    s <- sequences[i]
    if (!is.null(prob))
      prob <- prob[i, ]
    start <- sample(seq_len(nchar(s) - width + 1), 1, prob=prob)
    start
  })
  ends <- starts + width - 1
  stopifnot(any(ends <= nchar(sequences[1])))
  seqs <- substr(sequences, starts, ends)
  DataFrame(starts=starts, sequences=seqs)
}

scoreSegment <- function(segment, pfm, background) {
  stopifnot(nchar(segment) == ncol(pfm))
  segment.letters <- strsplit(segment, "")[[1]]
  p.background <- background/sum(background)
  pwm <- prop.table(pfm, 2)
  q <- prod(diag(pwm[segment.letters, seq_len(nchar(segment))]))
  p <- prod(p.background[segment.letters])
  q/p
}

samplingStep <- function(sequence, pfm, background, width=8) {
  positions <- seq_len(nchar(sequence[1]) - width + 1)
  segments <- sapply(positions, function(i) {
    substr(sequence, i, i + width - 1)
  })
  scores <- sapply(segments, function(s) {
    scoreSegment(s, pfm, background)
  })
  d <- data.frame(positions, segments, scores,
                  norm.score=scores/sum(scores),
                  row.names=NULL)
  d
}

s <- read.DNAStringSet("data/motif.fasta")[1:20]
 

## real motif: GGTTATAT
width <- 8
max.iter <- 60
iter <- 0
converge <- FALSE
s.pfm <- initPFM()
s.background <- initBackground()
s.background <- updateBackground(s, s.background)

## assumes equal sequence lengths
possible.starts <- seq_len(nchar(s[1]) - width + 1)

position.probs <- matrix(0, length(s), length(possible.starts),
                         dimnames=list(seq_along(s), possible.starts))

## the big loop
scores <- NULL
while (!converge && iter <= max.iter) {
  ## Sample positions/motifs in each sequence
  if (is.null(scores))
    a <- samplePositions(s)
  else
    a <- samplePositions(s, prob=position.probs)
  
  ## Add each random width sequence (starting motifs) to PWM
  for (i in seq_along(a$sequences)) {
    s.pfm <- updatePFM(a$sequences[i], s.pfm)
  }

  ## For all sequences, shift its motif to background and remove from
  ## PWM
  k <- sample(seq_along(a$sequence))
  for (i in k) {
    s.pfm.alt <- removeFromPWM(a$sequence[k], s.pfm)
    s.background.alt <- updateBackground(DNAStringSet(a$sequence[k]), s.background)
    scores <- samplingStep(s[k], s.pfm.alt, s.background.alt)
    position.probs[k, ] <- scores$norm.score
  }
  iter <- iter + 1
  
  message(sprintf("On iteration %d.", iter))
}

