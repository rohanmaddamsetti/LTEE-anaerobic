## aerobic-anaerobic-metagenomics.R by Rohan Maddamsetti.

## IMPORTANT TODO: numbers of aerobic and anaerobic genes don't exactly match those in
## measureTargetSize.py. Debug this before publication.

## count the cumulative number of stars over time, and plot.

library(tidyverse)

##########################################################################
## get the lengths of all genes in REL606. Do by running:
## python printEcoliIDs.py -i ../data/REL606.7.gbk > ../results/REL606_IDs.csv.
REL606.genes <- read.csv('../results/REL606_IDs.csv',as.is=TRUE) %>%
mutate(length=strtoi(length))

## get anaerobic-specific and aerobic-specific genes
## (written out by ArcAnalysisScript.R)
aerobic.genes <- read.csv('../results/aerobic-specific-genes.csv', as.is=TRUE)
anaerobic.genes <- read.csv('../results/anaerobic-specific-genes.csv', as.is=TRUE)

#### get the names of the randomly chosen gene sets (from measureTargetSize.py.)
## of course these are NOT aerobic or anaerobic genes. Just that cardinality is preserved.
random.aerobic.genes <- read.csv('../results/random-aerobic-set.csv',as.is=TRUE)
random.anaerobic.genes <- read.csv('../results/random-anaerobic-set.csv',as.is=TRUE)

## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv

mutation.data <- read.csv('../results/LTEE-metagenome-mutations.csv',
                          header=TRUE,as.is=TRUE) %>%
mutate(anaerobic=(Gene %in% anaerobic.genes$gene)) %>%
mutate(aerobic=(Gene %in% aerobic.genes$gene)) %>%
mutate(random.anaerobic=(Gene %in% random.anaerobic.genes$gene)) %>%
mutate(random.aerobic=(Gene %in% random.aerobic.genes$gene)) %>%
mutate(Generation=t0/10000)

##########################################################################
## look at mutations in the metagenomics dataset that occur
## in arcA binding sites.
get.arcA.motif.muts <- function(arcA.binding.sites,mutation.data) {

  filter.arcA.motif.muts <- function(arcA.motif) {
    arcA.motif.muts <- mutation.data %>%
    filter(Position >= arcA.motif$REL606_start) %>%
    filter(Position <= arcA.motif$REL606_end)
    return(arcA.motif.muts)
  }
  
  ## split arcA.binding.sites on each row, then filter mutation.data on those
  ## coordinates as output, then join those dataframes.
  arcA.motif.mutations <- arcA.binding.sites %>%
  droplevels() %>% ## don't map/reduce empty subsets
  split(.$K12_arcA_motif) %>%
  map_dfr(.f=filter.arcA.motif.muts)

  return(arcA.motif.mutations)
}

arcA.binding.sites <- read.csv('../results/K12_arcA_motifs_in_REL606.csv',as.is=TRUE)
arcA.motif.mutations <- get.arcA.motif.muts(arcA.binding.sites,mutation.data)
arcA.motif.fixations <- filter(arcA.motif.mutations,fixation==1)
## write out a table for Nkrumah.
write.csv(file='../results/LTEE_metagenome_arcA_motif_mutations.csv',arcA.motif.mutations)
##########################################################################
## calculate the probability of a given configuration of parallel mutations
## in the metagenomes across populations.
## This seems to be a neat test for positive selection!
## Seems to give the same answer as results in Tenaillon Nature paper,
## and unfortunately not effective for historical contingency per se.

## Keep this code for now, in case I get some idea that builds on this.
gene.multinom.probability <- function (mutation.data,gene,logp=TRUE) {
  population.probs <- mutation.data %>% group_by(Population) %>%
  summarize(total.muts=n()) %>% mutate(prob=total.muts/sum(total.muts))

  gene.data <- filter(mutation.data,Gene==gene,Annotation=='missense')

  ## if no mutations, return NA.
  if (nrow(gene.data) == 0) return(NA)
  
  gene.data.counts <- gene.data %>% group_by(Population) %>%
  summarize(muts=n())

  ## hack to get lengths of vector right for dmultinom.
  full.pop.column <- data.frame(Population=population.probs$Population)
  stat.df <- full_join(full.pop.column,gene.data.counts)
  stat.df[is.na(stat.df)] <- 0

  dmultinom(stat.df$muts, size = sum(stat.df$muts), prob=population.probs$prob, log = logp)
}

## look at a couple genes for a sanity check.
gene.multinom.probability(mutation.data,'arcB',FALSE)
gene.multinom.probability(mutation.data,'hslU',FALSE)

## draw every gene in the genome. calculate log probability and rank.
genes.vec <- unique(mutation.data$Gene)
prob.vec <- sapply(genes.vec,function(gene) {gene.multinom.probability(mutation.data,gene)})

multinom.result <- data.frame(Gene=genes.vec,Prob=prob.vec) %>% drop_na() %>%
arrange(Prob) %>% mutate(index=1:n()) %>%
mutate(anaerobic=(Gene %in% anaerobic.genes$gene)) %>%
mutate(aerobic=(Gene %in% aerobic.genes$gene))

## take a look at these distributions. The top genes are all under strong
## positive selection in the LTEE, as reported by Tenaillon et al.
multinom.plot <- ggplot(multinom.result,aes(x=index,y=-Prob,label=Gene,color)) + geom_point() + theme_classic()

multinom.plot2 <- ggplot(filter(multinom.result,anaerobic==TRUE),aes(x=index,y=-Prob,label=Gene,color)) + geom_point() + theme_classic()

multinom.plot3 <- ggplot(filter(multinom.result,aerobic==TRUE),aes(x=index,y=-Prob,label=Gene,color)) + geom_point() + theme_classic()

##########################################################################
## REALLY DIRTY HACK TO GET A QUICK ANSWER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## (but use this idea for refactoring: feed in as input here?)
##aerobic.genes <- random.aerobic.genes
##anaerobic.genes <- random.anaerobic.genes

####### END HACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## from measureIntergenicTargetSize.py:
## Length of intergenic regions: 487863
intergenic.length <- 487863

## from ArcAnalysisScript.R: ########
######################NOTE: DOUBLE CHECK CONSISTENT WITH measureTargetSize.py. output!!!!!
## anaerobic-specific gene length: 457593
## aerobic-specific gene length: 219286
anaerobic.gene.length <- 457593
aerobic.gene.length <- 219286

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

random.anaerobic.length <- filter(target.size.numbers,
                                  set=='random_anaerobic')$total_gene_length
random.anaerobic.synon.sites <- filter(target.size.numbers,
                                       set=='random_anaerobic')$synon_sites
random.anaerobic.nonsynon.sites <- filter(target.size.numbers,
                                          set=='random_anaerobic')$non_synon_sites


random.aerobic.length <- filter(target.size.numbers,
                                set=='random_aerobic')$total_gene_length
random.aerobic.synon.sites <- filter(target.size.numbers,
                                     set=='random_aerobic')$synon_sites
random.aerobic.nonsynon.sites <- filter(target.size.numbers,
                                        set=='random_aerobic')$non_synon_sites

total.length <- filter(target.size.numbers,set=='genome')$total_gene_length
total.synon.sites <- filter(target.size.numbers,set=='genome')$synon_sites
total.nonsynon.sites <- filter(target.size.numbers,set=='genome')$non_synon_sites

## DIRTY HACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##anaerobic.synon.sites <- random.anaerobic.synon.sites
##anaerobic.nonsynon.sites <- random.anaerobic.nonsynon.sites
##aerobic.synon.sites <- random.aerobic.synon.sites
##aerobic.nonsynon.sites <- random.aerobic.nonsynon.sites

##anaerobic.synon.sites <- random.anaerobic.length
##aerobic.synon.sites <- random.aerobic.length
##anaerobic.nonsynon.sites <- random.anaerobic.length
##aerobic.nonsynon.sites <- random.aerobic.length
##total.synon.sites <- total.length
##total.nonsynon.sites <- total.length

##anaerobic.gene.length <- random.anaerobic.length
##aerobic.gene.length <- random.aerobic.length

############# END HACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


intergenic.mutation.data <- filter(mutation.data,
                                   Gene=='intergenic')

dN.mutation.data <- filter(mutation.data,
                           Annotation=='missense')
dS.mutation.data <- filter(mutation.data,
                           Annotation=='synonymous')

## NOTE: rate of accumulation depends on the size of the gene set,
## when looking at the left and right tails.

## also compare genes in the tails of the multinomial hit distribution.
## left tail is strong selection. right tail fits the null well.

## ODD! way more mutations in top 100 or top 200 genes.
## but way fewer in the top 50! Take a look at 1-50, and 51-100 genes.

left.tail.genes <- multinom.result$Gene[1:50]
left.tail2.genes <- multinom.result$Gene[51:100]
##left.tail.genes <- multinom.result$Gene[1:100]
##left.tail.genes <- multinom.result$Gene[1:200]
left.tail.mutation.data <- filter(mutation.data,
                                  Gene %in% left.tail.genes)

left.tail2.mutation.data <- filter(mutation.data,
                                  Gene %in% left.tail2.genes)


right.tail.genes <- multinom.result$Gene[3564:3714]
##right.tail.genes <- multinom.result$Gene[3514:3714]
##right.tail.genes <- multinom.result$Gene[3614:3714]
right.tail.mutation.data <- filter(mutation.data,
                                   Gene %in% right.tail.genes)

## length of tails of multinomial distibution.
left.tail.length <- sum(filter(REL606.genes,gene %in% left.tail.genes)$length,na.rm=TRUE)
left.tail2.length <- sum(filter(REL606.genes,gene %in% left.tail2.genes)$length,na.rm=TRUE)
right.tail.length <- sum(filter(REL606.genes,gene %in% right.tail.genes)$length,na.rm=TRUE)


aerobic.mutation.data <- filter(mutation.data,aerobic==TRUE)
aerobic.dN.mutation.data <- filter(dN.mutation.data,aerobic==TRUE)
aerobic.dS.mutation.data <- filter(dS.mutation.data,aerobic==TRUE)

anaerobic.mutation.data <- filter(mutation.data,anaerobic==TRUE)
anaerobic.dN.mutation.data <- filter(dN.mutation.data,anaerobic==TRUE)
anaerobic.dS.mutation.data <- filter(dS.mutation.data,anaerobic==TRUE)


high.freq.dN.mutation.data <- filter(dN.mutation.data,final_frequency>0.8)
aerobic.high.freq.dN.mutation.data <- filter(high.freq.dN.mutation.data,aerobic==TRUE)
anaerobic.high.freq.dN.mutation.data <- filter(high.freq.dN.mutation.data,
                                               anaerobic==TRUE)

high.freq.dS.mutation.data <- filter(dS.mutation.data,final_frequency>0.8)
aerobic.high.freq.dS.mutation.data <- filter(high.freq.dS.mutation.data,aerobic==TRUE)
anaerobic.high.freq.dS.mutation.data <- filter(high.freq.dS.mutation.data,
                                               anaerobic==TRUE)

extinct.mutation.data <- filter(mutation.data,final_frequency==0)
extinct.dN.mutation.data <- filter(mutation.data,
                                   Annotation=='missense',
                                   final_frequency==0)
extinct.dS.mutation.data <- filter(mutation.data,
                                   Annotation=='synonymous',
                                   final_frequency==0)
extinct.aerobic.mutation.data <- filter(extinct.mutation.data,aerobic==TRUE)
extinct.anaerobic.mutation.data <- filter(extinct.mutation.data,anaerobic==TRUE)
extinct.aerobic.dN.mutation.data <- filter(extinct.dN.mutation.data,aerobic==TRUE)
extinct.anaerobic.dN.mutation.data <- filter(extinct.dN.mutation.data,anaerobic==TRUE)
extinct.aerobic.dS.mutation.data <- filter(extinct.dS.mutation.data,aerobic==TRUE)
extinct.anaerobic.dS.mutation.data <- filter(extinct.dS.mutation.data,anaerobic==TRUE)


high.freq.aerobic.muts <- filter(aerobic.mutation.data,final_frequency>0.8)
high.freq.anaerobic.muts <- filter(anaerobic.mutation.data,final_frequency>0.8)

c.high.freq.aerobic.mutations <- calc.cumulative.muts(high.freq.aerobic.muts,
                                            aerobic.gene.length)
c.high.freq.anaerobic.mutations <- calc.cumulative.muts(high.freq.anaerobic.muts,
                                                        anaerobic.gene.length)

c.anaerobic.mutations.no.dS <- calc.cumulative.muts(anaerobic.mutation.data.no.dS,
                                                    anaerobic.gene.length)
c.aerobic.mutations.no.dS <- calc.cumulative.muts(aerobic.mutation.data.no.dS,
                                                  aerobic.gene.length)


c.left.tail.mutations <- calc.cumulative.muts(left.tail.mutation.data,left.tail.length)
c.left.tail2.mutations <- calc.cumulative.muts(left.tail2.mutation.data,left.tail2.length)
c.right.tail.mutations <- calc.cumulative.muts(right.tail.mutation.data,right.tail.length)
c.intergenic.mutations <- calc.cumulative.muts(intergenic.mutation.data,intergenic.length)

c.mutations <- calc.cumulative.muts(mutation.data,total.length)
c.dN.mutations <- calc.cumulative.muts(dN.mutation.data,total.nonsynon.sites)
c.dS.mutations <- calc.cumulative.muts(dS.mutation.data,total.synon.sites)

c.high.freq.dS.mutations <- calc.cumulative.muts(high.freq.dS.mutation.data,
                                                 total.synon.sites) 

c.aerobic.mutations <- calc.cumulative.muts(aerobic.mutation.data,
                                            aerobic.gene.length)
c.anaerobic.mutations <- calc.cumulative.muts(anaerobic.mutation.data,
                                              anaerobic.gene.length)

c.aerobic.high.freq.dS.mutations <- calc.cumulative.muts(aerobic.high.freq.dS.mutation.data,
                                                         aerobic.synon.sites) 
c.anaerobic.high.freq.dS.mutations <- calc.cumulative.muts(anaerobic.high.freq.dS.mutation.data,
                                                           anaerobic.synon.sites) 
c.high.freq.dN.mutations <- calc.cumulative.muts(high.freq.dN.mutation.data,
                                                 total.nonsynon.sites) 
c.aerobic.high.freq.dN.mutations <- calc.cumulative.muts(aerobic.high.freq.dN.mutation.data,
                                                         aerobic.nonsynon.sites) 
c.anaerobic.high.freq.dN.mutations <- calc.cumulative.muts(anaerobic.high.freq.dN.mutation.data,
                                                           anaerobic.nonsynon.sites) 

c.extinct.mutations <- calc.cumulative.muts(extinct.mutation.data,
                                                             total.length)
c.extinct.dN.mutations <- calc.cumulative.muts(extinct.dN.mutation.data,
                                                             total.length)
c.extinct.dS.mutations <- calc.cumulative.muts(extinct.dS.mutation.data,
                                                             total.length)

c.extinct.aerobic.mutations <- calc.cumulative.muts(extinct.aerobic.mutation.data,
                                                             aerobic.gene.length)
c.extinct.anaerobic.mutations <- calc.cumulative.muts(extinct.anaerobic.mutation.data,
                                                             anaerobic.gene.length)
c.extinct.aerobic.dN.mutations <- calc.cumulative.muts(extinct.aerobic.dN.mutation.data,
                                                             aerobic.nonsynon.sites)
c.extinct.anaerobic.dN.mutations <- calc.cumulative.muts(extinct.anaerobic.dN.mutation.data,
                                                             anaerobic.nonsynon.sites)
c.extinct.aerobic.dS.mutations <- calc.cumulative.muts(extinct.aerobic.dS.mutation.data,
                                                             aerobic.synon.sites)
c.extinct.anaerobic.dS.mutations <- calc.cumulative.muts(extinct.anaerobic.dS.mutation.data,
                                                             anaerobic.synon.sites)

c.aerobic.dN.mutations <- calc.cumulative.muts(aerobic.dN.mutation.data,
                                               aerobic.nonsynon.sites)
c.aerobic.dS.mutations <- calc.cumulative.muts(aerobic.dS.mutation.data,
                                               aerobic.synon.sites)

c.anaerobic.dN.mutations <- calc.cumulative.muts(anaerobic.dN.mutation.data,
                                                 anaerobic.nonsynon.sites)
c.anaerobic.dS.mutations <- calc.cumulative.muts(anaerobic.dS.mutation.data,
                                                 anaerobic.synon.sites)

plot.aerobic.vs.anaerobic.muts <- function(aerobic.data,anaerobic.data,all.data) {
  ggplot(aerobic.data,aes(x=Generation,y=log(normalized.cs))) +
  theme_classic() +
  geom_point(size=0.2) +
  facet_wrap(.~Population,scales='fixed') +
  geom_point(data=anaerobic.data,aes(x=Generation,y=log(normalized.cs)),color='red',size=0.2) +
  geom_point(data=all.data,aes(x=Generation,y=log(normalized.cs)),color='grey',size=0.2) +
  ylab('Cumulative number of mutations, normalized by target size') +
  xlab('Generations (x 10,000)')
}

plot.aerobic.vs.anaerobic.muts2 <- function(aerobic.data,
                                            anaerobic.data,
                                            all.data,
                                            intergenic.data,
                                            left.tail.data,
                                            right.tail.data) {
  ggplot(aerobic.data,aes(x=Generation,y=log(normalized.cs))) +
  theme_classic() +
  geom_point(size=0.2) +
  facet_wrap(.~Population,scales='fixed') +
  geom_point(data=anaerobic.data,aes(x=Generation,y=log(normalized.cs)),color='red',size=0.2) +
  geom_point(data=all.data,aes(x=Generation,y=log(normalized.cs)),color='grey',size=0.2) +
  geom_point(data=intergenic.data,aes(x=Generation,y=log(normalized.cs)),color='orange',size=0.2) +
  geom_point(data=left.tail.data,aes(x=Generation,y=log(normalized.cs)),color='purple',size=0.2) +
  geom_point(data=right.tail.data,aes(x=Generation,y=log(normalized.cs)),color='blue',size=0.2) +
  ylab('Cumulative number of mutations, normalized by target size') +
  xlab('Generations (x 10,000)')
}

plot.aerobic.vs.anaerobic.muts3 <- function(aerobic.data,
                                            anaerobic.data,
                                            all.data,
                                            intergenic.data,
                                            left.tail.data,
                                            right.tail.data) {
  ggplot(aerobic.data,aes(x=Generation,y=normalized.cs)) +
  theme_classic() +
  geom_point(size=0.2) +
  facet_wrap(.~Population,scales='fixed') +
  geom_point(data=anaerobic.data,aes(x=Generation,y=normalized.cs),color='red',size=0.2) +
  geom_point(data=all.data,aes(x=Generation,y=normalized.cs),color='grey',size=0.2) +
  geom_point(data=intergenic.data,aes(x=Generation,y=normalized.cs),color='orange',size=0.2) +
  geom_point(data=left.tail.data,aes(x=Generation,y=normalized.cs),color='purple',size=0.2) +
  geom_point(data=right.tail.data,aes(x=Generation,y=normalized.cs),color='blue',size=0.2) +
  ylab('Cumulative number of mutations, normalized by target size') +
  xlab('Generations (x 10,000)')
}


## plot cumulative plot of anaerobic dN and aerobic dN on top of each other.
c.dN.plot <- plot.aerobic.vs.anaerobic.muts(c.aerobic.dN.mutations, c.anaerobic.dN.mutations,c.dN.mutations)
ggsave(c.dN.plot,filename='../results/figures/just-dN.pdf')

c.dS.plot <- plot.aerobic.vs.anaerobic.muts(c.aerobic.dS.mutations, c.anaerobic.dS.mutations,c.dS.mutations)
ggsave(c.dS.plot,filename='../results/figures/just-dS.pdf')

## This is weird! Check out if pattern holds if I split up by final frequency at 60K generations.
## The curves do seem to change depending on final frequency.

c.high.freq.dS.plot <- plot.aerobic.vs.anaerobic.muts(c.aerobic.high.freq.dS.mutations,c.anaerobic.high.freq.dS.mutations)
ggsave(c.high.freq.dS.plot,filename='../results/figures/highfreq-dS.pdf')

c.high.freq.dN.plot <- plot.aerobic.vs.anaerobic.muts(c.aerobic.high.freq.dN.mutations,c.anaerobic.high.freq.dN.mutations)
ggsave(c.high.freq.dN.plot,filename='../results/figures/highfreq-dN.pdf')

### This is interesting! Plot trajectory of occurrence of dN and dS in the whole population!
c.total.dN.dS.plot <- ggplot(c.dN.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point(color='purple') +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.dS.mutations,aes(x=Generation,y=normalized.cs,color='dS'),color='green')
ggsave(c.total.dN.dS.plot,filename='../results/figures/alldNdS.pdf')

## based on redoing the binomial test on genomes, it seems the difference in result
## is NOT driven by target size. rather, it is whether or not all kinds of mutations
## are included, and not just restricting to point mutations.

c.no.dS.plot <- plot.aerobic.vs.anaerobic.muts(c.aerobic.mutations.no.dS, c.anaerobic.mutations.no.dS)
ggsave(c.no.dS.plot,filename='../results/figures/all-mutations-but-dS.pdf')

##c.all.mut.plot <- plot.aerobic.vs.anaerobic.muts(c.aerobic.mutations, c.anaerobic.mutations,c.mutations)

c.all.mut.plot2 <- plot.aerobic.vs.anaerobic.muts2(c.aerobic.mutations, c.anaerobic.mutations,c.mutations,c.intergenic.mutations,c.left.tail.mutations,c.right.tail.mutations)
ggsave(c.all.mut.plot2,filename='../results/figures/all-mutations2.pdf')

c.all.mut.plot3 <- plot.aerobic.vs.anaerobic.muts3(c.aerobic.mutations, c.anaerobic.mutations,c.mutations,c.intergenic.mutations,c.left.tail.mutations,c.right.tail.mutations)
ggsave(c.all.mut.plot3,filename='../results/figures/all-mutations3.pdf')

## does this plot depend on final frequency?
c.high.freq.all.mut.plot <- plot.aerobic.vs.anaerobic.muts(c.high.freq.aerobic.mutations, c.high.freq.anaerobic.mutations)
ggsave(c.high.freq.all.mut.plot,filename='../results/figures/highfreq-all-muts.pdf')

### let's plot just extinct mutations. Any insights from this?
c.extinct.mut.plot <- plot.aerobic.vs.anaerobic.muts(c.extinct.aerobic.mutations, c.extinct.anaerobic.mutations,c.extinct.mutations)
ggsave(c.extinct.mut.plot,filename='../results/figures/extinct-all-muts.pdf')

c.extinct.dN.plot <- plot.aerobic.vs.anaerobic.muts(c.extinct.aerobic.dN.mutations,c.extinct.anaerobic.dN.mutations,c.extinct.dN.mutations)
ggsave(c.extinct.dN.plot,filename='../results/figures/extinct-dN.pdf')

c.extinct.dS.plot <- plot.aerobic.vs.anaerobic.muts(c.extinct.aerobic.dS.mutations,c.extinct.anaerobic.dS.mutations,c.extinct.dS.mutations)
ggsave(c.extinct.dS.plot,filename='../results/figures/extinct-dS.pdf')
