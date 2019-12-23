## aerobic-anaerobic-metagenomics.R by Rohan Maddamsetti.

## IMPORTANT TODO: numbers of aerobic and anaerobic genes don't exactly match those in
## measureTargetSize.py. Debug this before publication.

## Basic premise.
## count the cumulative number of stars over time, and plot.
## examine different kinds of mutations and genes.

## MAKE SURE THAT ZEROS IN MUTATION COUNTS ARE NOT THROWN OUT BY DPLYR!

library(tidyverse)
##########################################################################

## get the lengths of all genes in REL606.
## This excludes genes in repetitive regions of the genome.
## See Section 4.3.1 "Removing mutations in repetitive regions of the genome"
## in Ben Good's LTEE metagenomics paper for more details.
## This filtering is done in my python script printEcoliIDs.py.
##Do by running:
##python printEcoliIDs.py -i ../data/REL606.7.gbk > ../results/REL606_IDs.csv.
REL606.genes <- read.csv('../results/REL606_IDs.csv',as.is=TRUE) %>%
mutate(gene_length=strtoi(gene_length))

## get anaerobic-specific and aerobic-specific genes
## (written out by ArcAnalysisScript.R)
aerobic.genes <- read.csv('../results/aerobic-specific-genes.csv', as.is=TRUE)
anaerobic.genes <- read.csv('../results/anaerobic-specific-genes.csv', as.is=TRUE)

## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv

mutation.data <- read.csv('../results/LTEE-metagenome-mutations.csv',
                          header=TRUE,as.is=TRUE) %>%
    mutate(Generation=t0/10000) %>%
    mutate(anaerobic=(Gene %in% anaerobic.genes$gene)) %>%
    mutate(aerobic=(Gene %in% aerobic.genes$gene))

gene.mutation.data <- inner_join(mutation.data,REL606.genes)

## It turns out that some gene names map to multiple genes!!!
duplicate.genes <- gene.mutation.data %>%
    filter(Gene!='intergenic') %>%
    group_by(Gene) %>%
    summarize(checkme=length(unique(gene_length))) %>%
    filter(checkme>1)

## filter those duplicates.
gene.mutation.data <- filter(gene.mutation.data, !(Gene %in% duplicate.genes$Gene))

####### Constants. REVAMP THIS CODE!
## from measureIntergenicTargetSize.py:
## Length of intergenic regions: 487863
intergenic.length <- 487863

## from aerobic-anaerobic-genomics.R: ########
######################NOTE: DOUBLE CHECK CONSISTENT WITH measureTargetSize.py. output!!!!!

anaerobic.gene.length <- 457593
aerobic.gene.length <- 219286
##########################################################################

## Normalization constant calculations.
## IMPORTANT-- THIS STEP IS REALLY IMPORTANT.
## NOT CLEAR HOW TO NORMALIZE IN ORDER TO COMPARE DIFFERENT CLASSES OF MUTATIONS
## "APPLES TO APPLES".


## TODO: check difference between normalizing by gene length and normalizing
## by synonymous/nonsynonymous opportunities. The end result should be very similar.
## Then consider refactoring to just use REL606_IDs.csv to get gene lengths.

## numbers gotten by running measureTargetSize.py.
## Use these to normalize cumulative mutations over time.
target.size.numbers <- read.csv('../results/target_size.csv',header=TRUE,as.is=TRUE)

anaerobic.length <- filter(target.size.numbers,set=='anaerobic')$total_gene_length
anaerobic.synon.sites <- filter(target.size.numbers,set=='anaerobic')$synon_sites
anaerobic.nonsynon.sites <- filter(target.size.numbers,set=='anaerobic')$non_synon_sites

aerobic.length <- filter(target.size.numbers,set=='aerobic')$total_gene_length
aerobic.synon.sites <- filter(target.size.numbers,set=='aerobic')$synon_sites
aerobic.nonsynon.sites <- filter(target.size.numbers,set=='aerobic')$non_synon_sites


total.length <- filter(target.size.numbers,set=='genome')$total_gene_length
total.synon.sites <- filter(target.size.numbers,set=='genome')$synon_sites
total.nonsynon.sites <- filter(target.size.numbers,set=='genome')$non_synon_sites

########################################
## look at accumulation of stars over time.
## in other words, look at the rates at which the mutations occur over time.
## plot cumulative sum of anaerobic and aerobic dS and dN in each population.
cumsum.per.pop.helper.func <- function(df) {
  df %>%
  arrange(t0) %>%
  group_by(Population,Generation) %>%
  summarize(count=n()) %>%
  mutate(cs=cumsum(count))
}

calc.cumulative.muts <- function(data, normalization.constant) {
  data %>%
  split(.$Population) %>%
  map_dfr(.f=cumsum.per.pop.helper.func) %>%
  mutate(normalized.cs=cs/normalization.constant)
}
#########################################################################
## calculate cumulative numbers of mutations in each category.
## Then, make plots to see if any interesting patterns emerge.

plot.cumulative.muts <- function(mut.data,logscale=TRUE, my.color="black") {
    if (logscale) {
        p <- ggplot(mut.data,aes(x=Generation,y=log(normalized.cs)))
    } else {
        p <- ggplot(mut.data,aes(x=Generation,y=normalized.cs))
    }
    p <- p +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        geom_step(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='fixed') +
        ylab('Cumulative number of mutations, normalized by target size') +
        xlab('Generations (x 10,000)')
    return(p)
}

## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
add.cumulative.mut.layer <- function(p, layer.df, my.color, logscale=TRUE) {
    if (logscale) {
        p <- p +
            geom_point(data=layer.df, aes(x=Generation,y=log(normalized.cs)), color=my.color, size=0.2) +
            geom_step(data=layer.df, size=0.2, color=my.color)
        } else {
            p <- p +
                geom_point(data=layer.df, aes(x=Generation,y=normalized.cs), color=my.color, size=0.2) +
                geom_step(data=layer.df, size=0.2, color=my.color)
        }
    return(p)
}
#########################################################################
## examine dS over the genome.
gene.dS.mutation.data <- gene.mutation.data %>%
filter(Annotation=='synonymous')

## examing dN over the genome.
gene.dN.mutation.data <- gene.mutation.data %>%
filter(Annotation=='missense')

## let's look at nonsense mutations.
gene.nonsense.mutation.data <- gene.mutation.data %>%
filter(Annotation=='nonsense')

## let's look at all mutations.
c.mutations <- calc.cumulative.muts(mutation.data,total.length)

##c.dN.mutations <- calc.cumulative.muts(dN.mutation.data,total.nonsynon.sites)
dN.normalization.const <- total.length * total.nonsynon.sites/(total.synon.sites+total.nonsynon.sites)
c.dN.mutations <- calc.cumulative.muts(gene.dN.mutation.data,dN.normalization.const)

##c.dS.mutations <- calc.cumulative.muts(dS.mutation.data,total.synon.sites)
dS.normalization.const <- total.length * total.synon.sites/(total.synon.sites+total.nonsynon.sites)
c.dS.mutations <- calc.cumulative.muts(gene.dS.mutation.data, dS.normalization.const)


aerobic.mutation.data <- filter(mutation.data,aerobic==TRUE)
aerobic.dN.mutation.data <- filter(aerobic.mutation.data,
                                   Annotation=='missense')
aerobic.dS.mutation.data <- filter(aerobic.mutation.data,
                                   Annotation=='synonymous')
aerobic.sv.mutation.data <- filter(aerobic.mutation.data,
                                   Annotation=='sv')
aerobic.indel.mutation.data <- filter(aerobic.mutation.data,
                                   Annotation=='indel')

anaerobic.mutation.data <- filter(mutation.data,anaerobic==TRUE)
anaerobic.dN.mutation.data <- filter(anaerobic.mutation.data,
                                     Annotation=='missense')
anaerobic.dS.mutation.data <- filter(anaerobic.mutation.data,
                                     Annotation=='synonymous')
anaerobic.sv.mutation.data <- filter(anaerobic.mutation.data,
                                   Annotation=='sv')
anaerobic.indel.mutation.data <- filter(anaerobic.mutation.data,
                                   Annotation=='indel')

c.aerobic.mutations <- calc.cumulative.muts(aerobic.mutation.data,
                                            aerobic.gene.length)
c.anaerobic.mutations <- calc.cumulative.muts(anaerobic.mutation.data,
                                              anaerobic.gene.length)

c.aerobic.sv.mutations <- calc.cumulative.muts(aerobic.sv.mutation.data,
                                            aerobic.gene.length)
c.anaerobic.sv.mutations <- calc.cumulative.muts(anaerobic.sv.mutation.data,
                                              anaerobic.gene.length)

c.aerobic.indel.mutations <- calc.cumulative.muts(aerobic.indel.mutation.data,aerobic.gene.length)
c.anaerobic.indel.mutations <- calc.cumulative.muts(anaerobic.indel.mutation.data,anaerobic.gene.length)

c.aerobic.dN.mutations <- calc.cumulative.muts(aerobic.dN.mutation.data,
                                               aerobic.nonsynon.sites)
c.aerobic.dS.mutations <- calc.cumulative.muts(aerobic.dS.mutation.data,
                                               aerobic.synon.sites)

c.anaerobic.dN.mutations <- calc.cumulative.muts(anaerobic.dN.mutation.data,
                                                 anaerobic.nonsynon.sites)
c.anaerobic.dS.mutations <- calc.cumulative.muts(anaerobic.dS.mutation.data,
                                                 anaerobic.synon.sites)

#############################################################################

## Figures.

## plot all classes of mutations in aerobic versus anaerobic genes.
c.aerobic.vs.anaerobic.plot1 <- plot.cumulative.muts(c.aerobic.mutations,my.color='black',logscale=TRUE) %>%
    add.cumulative.mut.layer(c.anaerobic.mutations, my.color='red',logscale=TRUE) %>%
    add.cumulative.mut.layer(c.mutations, my.color='grey',logscale=TRUE)

c.aerobic.vs.anaerobic.plot1


## all mutations in gray
sv.indel.nonsense.mutation.data <- mutation.data %>%
    filter(Annotation %in% c("sv", "indel", "nonsense"))
c.sv.indel.nonsense.mutations <- calc.cumulative.muts(sv.indel.nonsense.mutation.data,total.length)

## aerobic in black
aerobic.sv.indel.nonsense.muts <- mutation.data %>%
    filter(aerobic==TRUE) %>%
    filter(Annotation %in% c("sv", "indel", "nonsense"))
c.aerobic.sv.indel.nonsense.mutations <- calc.cumulative.muts(aerobic.sv.indel.nonsense.muts,aerobic.gene.length)

## anaerobic in red
anaerobic.sv.indel.nonsense.muts <- mutation.data %>%
    filter(anaerobic==TRUE) %>%
    filter(Annotation %in% c("sv", "indel", "nonsense"))
c.anaerobic.sv.indel.nonsense.mutations <- calc.cumulative.muts(anaerobic.sv.indel.nonsense.muts,anaerobic.gene.length)

aerobic.anaerobic.purifying.plot1 <- plot.cumulative.muts(
    c.aerobic.sv.indel.nonsense.mutations, logscale=TRUE) %>%
    add.cumulative.mut.layer(c.anaerobic.sv.indel.nonsense.mutations, my.color="red",
                             logscale=TRUE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsense.mutations, my.color="grey",
                             logscale=TRUE)
aerobic.anaerobic.purifying.plot1

aerobic.anaerobic.purifying.plot2 <- plot.cumulative.muts(
    c.aerobic.sv.indel.nonsense.mutations, logscale=FALSE) %>%
    add.cumulative.mut.layer(c.anaerobic.sv.indel.nonsense.mutations, my.color="red",
                             logscale=FALSE) %>%
    add.cumulative.mut.layer(c.sv.indel.nonsense.mutations, my.color="grey",
                             logscale=FALSE)
aerobic.anaerobic.purifying.plot2

## plot transposons and indels as separate layers.
c.sv.indels.plot <- plot.cumulative.muts(c.aerobic.sv.mutations,my.color='black',logscale=TRUE) %>%
    add.cumulative.mut.layer(c.aerobic.indel.mutations, my.color='brown',logscale=TRUE) %>%
    add.cumulative.mut.layer(c.anaerobic.sv.mutations, my.color='red',logscale=TRUE) %>%
    add.cumulative.mut.layer(c.anaerobic.indel.mutations, my.color='pink',logscale=TRUE)

c.sv.indels.plot

##########################################################################
## look at accumulation of stars over time for random subsets of genes.
## I want to plot the distribution of cumulative mutations over time for
## say, 1000 or 10000 random subsets of genes.

plot.random.subsets <- function(data, subset.size=300, N=1000,log=TRUE) {

  ## set up an empty plot then add random trajectories, one by one.
  my.plot <- ggplot(data) +
  theme_classic() +
  facet_wrap(.~Population,scales='fixed') +
  ylab('Cumulative number of mutations, normalized by target size') +
  xlab('Generations (x 10,000)')

  for (i in 1:N) {
    rando.genes <- sample(unique(data$Gene),subset.size)
    mut.subset <- filter(data,Gene %in% rando.genes)
    subset.length <- sum(mut.subset$length)
    c.mut.subset <- calc.cumulative.muts(mut.subset,subset.length)
    
    if (log) {
      my.plot <- my.plot +
      geom_point(data=c.mut.subset,aes(x=Generation,y=log(normalized.cs)), color='gray',size=0.2,alpha = 0.1)
    } else {
      my.plot <- my.plot +
      geom_point(data=c.mut.subset,aes(x=Generation,y=normalized.cs), color='gray',size=0.2,alpha = 0.1)
    }
  }
  return(my.plot)
}

log.rando.plot <- plot.random.subsets(full.mutation.data,log=TRUE)
ggsave(log.rando.plot,filename='../results/figures/log-rando-plot.png')
rando.plot <- plot.random.subsets(full.mutation.data,log=FALSE)
ggsave(rando.plot,filename='../results/figures/rando-plot.png')
