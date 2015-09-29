# PairDistances

PairDistances contains code for performing pairwise alignment of protein sequences. We avoid the use of an affine gap penalty by simply not inserting gaps. We find the most-likely alignment configuration and evolutionary time between each pair of sequences using a Markov Process approach (see Thorne - ``An Evolutionary model for Maximum Likelihood Alignment of DNA Sequences'').

The user simply needs to provide a FASTA file for the sequences in which they which to calculate the distances between and a results tab-delimited file will be generated. The format of the file follows that shown below:


| Seq1        | Seq2           | Time  | LogLikelihood | Configuration | S1Length | S2Length | S1Sequence | S2Sequence |
|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| P21538 | P21538 | 0.0 | -11.440250364624116 | 3.0 | 4 | 4 | KKRK | KKRK |
| P00546 | P21538 | 0.9582527115139087 | -32.41837403019155 | 5.0 | 7 | 4 | PQWRRKD | KKRK |
| P00546 | P00546 | 0.0 | -22.16493324085298 | 6.0 | 7 | 7 | PQWRRKD | PQWRRKD |
| P37261 | P21538 | 0.1473260089123268 | -17.801589031708815 | 3.0 | 4 | 4 | KKRP | KKRK |
| P37261 | P00546 | 1.7000000000000002 | -33.38128364660126 | 6.0 | 4 | 7 | KKRP | PQWRRKD |
| P37261 | P37261 | 0.0 | -11.603586264197588 | 3.0 | 4 | 4 | KKRP | KKRP |
| P25555 | P21538 | 1.1278615505309812 | -30.37755252356034 | 5.0 | 7 | 4 | PVRRRLS | KKRK |


There are a variety of input parameters that the user can provide to improve both the efficiency and accuracy of the pairwise alignment search. Details of which will be added here in the future, but in the mean time, look in the PairDistances.java file.
