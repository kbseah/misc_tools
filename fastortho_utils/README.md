# Working with FastOrtho output

## FastOrtho 
[FastOrtho](http://enews.patricbrc.org/fastortho/) is a reimplementation of the OrthoMCL pipeline for finding ortholog clusters using the [MCL](http://micans.org/mcl/) algorithm.

### Advantages:
 * Fast clustering algorithm, doesn't require SQL database (unlike OrthoMCL)
 * Avoids double-counting problem (e.g. early versions of LS-BSR)
 * Easy-to-use graphical interface for configuring run parameters 

### Disadvantages:
 * Blast process is not parallelized - is the rate-limiting bottleneck (but can be circumvented by pre-computing blast results and using that as input)

## Utility scripts

FastOrtho only finds ortholog clusters and output is in its own idiosyncratic format. These scripts are to help in parsing that output to find core/pan genomes etc.

### Extract single-copy orthologs

Extract ortholog clusters that are in single-copy in all genomes and align with Muscle (assumes that `muscle` binary is in path)

```bash
 perl fastortho2fasta.pl -e fastortho.end -n <num_taxa> -f fastortho.faa -c concatenated alignment.fasta
```

### Resampling for core/pan genome statistics

Generate accumulation curves of core/pan genome counts vs. number of genomes, with error bars derived from resampling (with replacement!). 

```bash
 perl resampling_FastOrtho.pl fastortho.end fastortho.opt > result.out
```

Results can be imported into R and plotted, e.g.:

```R
 d <- read.table("result.out",header=T)
 plot(d$num_genomes,d$pan,pch=".",ylim=c(0,max(d$pan)))
 points(d$num_genomes,d$core,pch=".",col="blue")
 points(d$num_genomes,d$singleton,pch=".",col="red")
 points(tapply(X=d$pan,INDEX=d$num_genomes,FUN=mean),type="l")
 points(tapply(X=d$core,INDEX=d$num_genomes,FUN=mean),type="l",col="blue")
 points(tapply(X=d$singleton,INDEX=d$num_genomes,FUN=mean),type="l",col="red")
 legend(x="topleft",col=c("black","red","blue"),legend=c("pan","singleton","core"),lty=1)
```

or with the `ggplot2` and `reshape2` packages:

```R
 d <- read.table("result.out",header=T)
 library(ggplot2)
 library(reshape2)
 d.melt <- melt(d,id.vars="num_genomes") # Reformat data to "long" fromat
 ggplot(d.melt,aes(num_genomes,value,col=variable)) + 
  geom_jitter(height=0,width=0.1,size=0.1) + # jitter in horizontal axis for legibility
  stat_summary(fun.y="mean",geom="line")
```
