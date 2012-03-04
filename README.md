# Toy Gibbs Sampling Motif Finder

This is a "toy" (meaning fun to play with, but *don't* use it as a
tool) Gibbs Sampling motif finder.

Semanatically speaking, it's actually a *site sampler*. Why? Motif
samplers allow for more than one occurrence of more than one or more
motifs; site samplers simply assume a motif occurs at only one site in
each sequence.

This was written in an afternoon; the code is a bit messy,
underdocumented, and ignores some ideas of the original paper (i.e. I
do away with psuedocounts entirely). That said, it works pretty well
(even though it's painfully slow).

## How does it work?

I'll explain this soon.

## How do I run it?

Just `source("R/motif.R")` in the project directory. There is *no
convergence checking* yet! It finds a motif in the toy data after
about 50 iterations.

## Why is it so slow?

Use `Rprof()` to see why. `samplingStep` is a culprit. I'll write a
more useful "toy" version that uses some `.Call`s to speed it up.

## License

GPL 3.0.


