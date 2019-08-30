## aerobic-anaerobic-metagenomics.R by Rohan Maddamsetti.

## IMPORTANT TODO: numbers of aerobic and anaerobic genes don't exactly match those in
## measureTargetSize.py. Debug this before publication.

## Basic premise.
## count the cumulative number of stars over time, and plot.
## examine different kinds of mutations and genes.

## can I control for total number of mutations? Include as a predictor or something?

library(tidyverse)

##########################################################################
## get the lengths of all genes in REL606. Do by running:
## python printEcoliIDs.py -i ../data/REL606.7.gbk > ../results/REL606_IDs.csv.
REL606.genes <- read.csv('../results/REL606_IDs.csv',as.is=TRUE) %>%
mutate(gene_length=strtoi(gene_length))

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

full.mutation.data <- inner_join(mutation.data,REL606.genes)
## It turns out that some gene names map to multiple genes!!!
duplicate.genes <- full.mutation.data %>%
filter(Gene!='intergenic') %>%
group_by(Gene) %>%
summarize(checkme=length(unique(gene_length))) %>%
filter(checkme>1)
## filter those duplicates.
full.mutation.data <- filter(full.mutation.data, !(Gene %in% duplicate.genes$Gene))


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
gene.multinom.probability <- function (mutation.data,gene,logp=TRUE,dS=FALSE) {
  population.probs <- mutation.data %>% group_by(Population) %>%
  summarize(total.muts=n()) %>% mutate(prob=total.muts/sum(total.muts))

  if (dS) {
      gene.data <- filter(mutation.data,Gene==gene,Annotation=='synonymous')
  } else {
      gene.data <- filter(mutation.data,Gene==gene,Annotation=='missense')
  }

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
dN.prob.vec <- sapply(genes.vec,function(gene) {gene.multinom.probability(mutation.data,gene)})
dS.prob.vec <- sapply(genes.vec,function(gene) {gene.multinom.probability(mutation.data,gene,dS=TRUE)})

multinom.result <- data.frame(Gene=genes.vec,dN.Prob=dN.prob.vec,dS.Prob=dS.prob.vec) %>% drop_na() %>%
arrange(dN.Prob) %>% mutate(index=1:n()) %>%
mutate(anaerobic=(Gene %in% anaerobic.genes$gene)) %>%
mutate(aerobic=(Gene %in% aerobic.genes$gene))

## how much is this result related to the number of mutations in each gene?
## Control for gene length, and see if there's a correlation between log(p)
## and density of mutations in each gene.
## there's a super strong relationship. p < 2.26e-16.

## just rank gene.mutation.density to see what it looks like!
## maybe this is already a good indication of selection.
gene.mutation.density <- full.mutation.data %>%
filter(Gene!='intergenic') %>%
group_by(Gene) %>%
summarize(density=n()/unique(gene_length)) %>%
arrange(desc(density))

multinom.result.and.density <- inner_join(gene.mutation.density,multinom.result)
cor.test(multinom.result.and.density$density,multinom.result.and.density$Prob)

########################################################################
## investigate dS across the genome in the metagenomics data.
## reuse code from my 2015 Mol. Biol. Evol. paper.

ks.analysis <- function (the.data) {
  ## For each set of data (all data, non-mutators, MMR mutators, mutT mutators)
  ## do the following: 1) make a uniform cdf on mutation rate per base.
  ## 2) make a thetaS cdf. 3) make an empirical cdf of mutations per gene.
  ## do K-S tests for goodness of fit of the empirical cdf with the cdfs for
  ## the uniform cdf and thetaS cdf hypotheses.

  hit.genes.df <- the.data %>%
  group_by(Gene,locus_tag,gene_length) %>%
  summarize(hits=n()) %>%
  ungroup() %>%
  arrange(desc(gene_length))
  
  hit.genes.length <- sum(hit.genes.df$gene_length)
  
  ## Calculate the empirical distribution of synonymous substitutions per gene.
  mutation.total <- sum(hit.genes.df$hits)
  empirical.cdf <- cumsum(hit.genes.df$hits)/mutation.total
  ## Null hypothesis: probability of a mutation per base is uniform.
  null.cdf <- cumsum(hit.genes.df$gene_length)/hit.genes.length)

  ## Do Kolmogorov-Smirnov tests for goodness of fit.
  print(ks.test(empirical.cdf, null.cdf, simulate.p.value=TRUE))

  results.to.plot <- data.frame(locus_tag=hit.genes.df$locus_tag, Gene=hit.genes.df$Gene, gene_length=hit.genes.df$gene_length, empirical=empirical.cdf,null=null.cdf)

  return(results.to.plot)
}

make.KS.Figure <- function(the.results.to.plot) {
## This function generates the first panel that I want
## but not with synonymous substitution rates in the LTEE.
  
  ## for plotting convienence, add an index to the data frame.
  the.results.to.plot$index <- 1:length(the.results.to.plot$gene)
  
  plot <- ggplot(the.results.to.plot, aes(x=index)) +
    geom_line(aes(y=empirical), colour="red") + 
    geom_line(aes(y=null), linetype=2) + 
    scale_x_continuous('Genes ranked by length',limits=c(0,2900)) +
    scale_y_continuous('Cumulative proportion of synonymous mutations',limits=c(0,1)) +
      theme_classic() + theme(axis.title=element_text(size=18),axis.text=element_text(size=12))
  plot
  
}

full.dS.mutation.data <- full.mutation.data %>%
filter(Gene!='intergenic') %>%
filter(Annotation=='synonymous')

## dS mutations occurring in metagenomics don't fall uniformly over the genome!
cumsum.dS.over.metagenome <- ks.analysis(full.dS.mutation.data)
make.KS.Figure(cumsum.dS.over.metagenome)


## cross check with dS in the mutator 50K genomes.
mutator.50K <- read.csv('../data/Gen50000_M.csv',header=TRUE,stringsAsFactors=FALSE) %>%
  select(population,strain,clone,mutator_status,type,start_position,end_position,gene_position,
         locus_tag,mutation_category,snp_type) %>% filter(clone=='A')

genome.colnames.to.metagenome.colnames <- function(genome.data) {
  renamed <- genome.data %>% rename(Population=population,
                                    Locus_tag = locus_tag)
  return(renamed)
}


mutator.50K.dS <- filter(mutator.50K,snp_type=='synonymous') %>%
genome.colnames.to.metagenome.colnames %>%
left_join(full.mutation.data)

cumsum.dS.over.genome <- ks.analysis(mutator.50K.dS)
## FASCINATING!!! Looks like dS is not neutral in the
## genomes either!!

make.KS.Figure(cumsum.dS.over.genome)

## IMPORTANT TODO: Double-check with the 50K and 40K genomes too, to see whether
## this result is due to a bug in my code somewhere!

mutator.40K <- read.csv('../data/Gen40000_M.csv',header=TRUE,as.is=TRUE) %>%
  select(population,strain,clone,mutator_status,type,start_position,end_position,gene_position,
         locus_tag,mutation_category,snp_type) %>% filter(clone=='A')


### do the same KS.test, but a different way to double-check.
second.test <- mutator.50K.dS %>% genome.colnames.to.metagenome.colnames %>%
group_by(Locus_tag) %>% summarize(dS.count=n()) %>%
inner_join(select(cumsum.dS.over.metagenome,Locus_tag=locus_tag,Gene=Gene,gene_length=gene_length)) %>%
arrange(gene_length) %>%
mutate(index=1:n()) %>%
mutate(dS.density=dS.count/gene_length) %>%
mutate(empirical=cumsum(dS.count)/sum(dS.count)) %>%
mutate(null=cumsum(gene_length)/sum(gene_length))



genome.dS.density.plot <- ggplot(second.test,aes(x=index,y=dS.density,label=Gene)) +
geom_point() +
theme_classic()

genome.dS.KS.plot <- ggplot(second.test,aes(x=index)) +
    geom_line(aes(y=empirical), colour="red") + 
    geom_line(aes(y=null), linetype=2) + 
    scale_x_continuous('Genes ranked by length',limits=c(0,2900)) +
    scale_y_continuous('Cumulative proportion of synonymous mutations',limits=c(0,1)) +
      theme_classic() + theme(axis.title=element_text(size=18),axis.text=element_text(size=12))
genome.dS.KS.plot
## This plot 


## This means that the dynamics of mutations under relaxed selection should be
## dragged by mutations under positive selection.

## Does the Ornstein-Uhlenbeck process require that effects of each mutation is uncorrelated?
## If mutations affecting anaerobic fitness are uncorrelated with mutations affecting
## aerobic fitness, does this 

## Genotype-Phenotype Map

## new null hypothesis: dynamics of ALL mutations are driven by positive selection.
## either hitchhiking or as a driver.


################################################################ EDITING HERE!

## are dS mutations concentrated in a few genes? Based on my 2015 paper, I expect a uniform
## distribution, after normalizing by gene length. Does the null for genomes work for
## the metagenome?

dS.mutation.density <- full.dS.mutation.data %>%
group_by(Gene,gene_length) %>%
summarize(dS.count=n()) %>%
ungroup() %>%
mutate(density=dS.count/gene_length)

dN.mutation.density <- full.mutation.data %>%
filter(Gene!='intergenic') %>%
filter(Annotation=='missense') %>%
group_by(Gene,gene_length) %>%
summarize(dN.count=n()) %>%
ungroup() %>%
mutate(density=dN.count/gene_length)

synon.density.plot <- ggplot(filter(dS.mutation.density,dS.count>3),aes(x=gene_length,y=density,label=Gene)) +
geom_point() +
theme_classic()

synon.count.plot <- ggplot(dS.mutation.density,aes(x=gene_length,y=dS.count,label=Gene,group=dS.count)) +
##geom_point() +
geom_violin() +
theme_classic()

dN.density.plot <- ggplot(filter(dN.mutation.density,dN.count>3),aes(x=gene_length,y=density,label=Gene)) +
geom_point() +
theme_classic()

test.model <- lm(data=dS.mutation.density,formula=dS.count~gene_length+0)

## don't know if this plot is informative or not...
## interesting that the really short gene ECB_01628 has two
## coexisting dS in Ara+6, but can't calculate significance post-hoc!
## genes in the left hand of the tail tend to be quite short.
synon.density.plot
synon.count.plot

##test <- full.mutation.data %>% filter(Gene=='ECB_01628')
test <- full.mutation.data %>% filter(Gene=='pmrD')

## write to file.
write.csv(multinom.result,file='../results/multinomial_test_for_selection.csv')

## take a look at these distributions. The top genes are all under strong
## positive selection in the LTEE, as reported by Tenaillon et al.
multinom.plot <- ggplot(multinom.result,aes(x=index,y=-Prob,label=Gene,color)) + geom_point() + theme_classic()

multinom.plot2 <- ggplot(filter(multinom.result,anaerobic==TRUE),aes(x=index,y=-Prob,label=Gene,color)) + geom_point() + theme_classic()

multinom.plot3 <- ggplot(filter(multinom.result,aerobic==TRUE),aes(x=index,y=-Prob,label=Gene,color)) + geom_point() + theme_classic()

## let's select the middle 2001 genes of the cumulative mutation graph, to plot
## an estimate of the median, rather than averaging over the whole genome,
## which could be skewed by the tails.

gene.mutation.density.midpoint <- round(length(gene.mutation.density$density)/2)
gmd.q2 <- gene.mutation.density.midpoint - 1000
gmd.q3 <- gene.mutation.density.midpoint + 1000

median.genes <- gene.mutation.density$Gene[gmd.q2:gmd.q3]
median.mutation.data <- filter(full.mutation.data,Gene %in% median.genes)
median.gene.length <- sum(median.mutation.data$gene_length)
c.median.mutations <- calc.cumulative.muts(median.mutation.data,median.gene.length)

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

intergenic.point.mutation.data <- filter(intergenic.mutation.data,Annotation=='noncoding')

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
c.intergenic.point.mutations <- calc.cumulative.muts(intergenic.point.mutation.data,intergenic.length)


c.mutations <- calc.cumulative.muts(mutation.data,total.length)
##c.dN.mutations <- calc.cumulative.muts(dN.mutation.data,total.nonsynon.sites)
dN.normalization.const <- total.length * total.nonsynon.sites/(total.synon.sites+total.nonsynon.sites)
c.dN.mutations <- calc.cumulative.muts(dN.mutation.data,dN.normalization.const)
##c.dS.mutations <- calc.cumulative.muts(dS.mutation.data,total.synon.sites)
dS.normalization.const <- total.length * total.synon.sites/(total.synon.sites+total.nonsynon.sites)
##c.dS.mutations <- calc.cumulative.muts(dS.mutation.data, total.synon.sites)
c.dS.mutations <- calc.cumulative.muts(dS.mutation.data, dS.normalization.const)

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
                                            dS.data
                                            ) {
  ggplot(aerobic.data,aes(x=Generation,y=log(normalized.cs))) +
  theme_classic() +
  geom_point(size=0.2) +
  facet_wrap(.~Population,scales='fixed') +
  geom_point(data=anaerobic.data,aes(x=Generation,y=log(normalized.cs)),color='red',size=0.2) +
  geom_point(data=all.data,aes(x=Generation,y=log(normalized.cs)),color='grey',size=0.2) +
  geom_point(data=intergenic.data,aes(x=Generation,y=log(normalized.cs)),color='orange',size=0.2) +
  geom_point(data=dS.data,aes(x=Generation,y=log(normalized.cs)),color='yellow',size=0.2) +
##  geom_point(data=left.tail.data,aes(x=Generation,y=log(normalized.cs)),color='purple',size=0.2) +
  ##geom_point(data=right.tail.data,aes(x=Generation,y=log(normalized.cs)),color='blue',size=0.2) +
  ylab('Cumulative number of mutations, normalized by target size') +
  xlab('Generations (x 10,000)')
}

plot.aerobic.vs.anaerobic.muts3 <- function(aerobic.data,
                                            anaerobic.data,
                                            all.data,
                                            intergenic.data,
                                            dS.data##,
                                            ##left.tail.data
                                            ) {
  ggplot(aerobic.data,aes(x=Generation,y=normalized.cs)) +
  theme_classic() +
  geom_point(size=0.2) +
  facet_wrap(.~Population,scales='fixed') +
  geom_point(data=anaerobic.data,aes(x=Generation,y=normalized.cs),color='red',size=0.2) +
  geom_point(data=all.data,aes(x=Generation,y=normalized.cs),color='grey',size=0.2) +
  geom_point(data=intergenic.data,aes(x=Generation,y=normalized.cs),color='orange',size=0.2) +
  geom_point(data=dS.data,aes(x=Generation,y=normalized.cs),color='yellow',size=0.2) +
##  geom_point(data=left.tail.data,aes(x=Generation,y=normalized.cs),color='purple',size=0.2) +
##  geom_point(data=right.tail.data,aes(x=Generation,y=normalized.cs),color='blue',size=0.2) +
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

c.all.mut.plot2 <- plot.aerobic.vs.anaerobic.muts2(c.aerobic.mutations, c.anaerobic.mutations,c.mutations,c.intergenic.mutations, c.dS.mutations)
ggsave(c.all.mut.plot2,filename='../results/figures/all-mutations2.pdf')

c.all.mut.plot3 <- plot.aerobic.vs.anaerobic.muts3(c.aerobic.mutations, c.anaerobic.mutations,c.mutations,c.intergenic.mutations, c.dS.mutations)
ggsave(c.all.mut.plot3,filename='../results/figures/all-mutations3.pdf')

######### Plot median instead of genome-wide average.
##c.all.mut.plot4 <- plot.aerobic.vs.anaerobic.muts2(c.aerobic.mutations, c.anaerobic.mutations,c.median.mutations,c.intergenic.mutations, c.dS.mutations)
c.all.mut.plot4 <- plot.aerobic.vs.anaerobic.muts2(c.aerobic.mutations, c.anaerobic.mutations,c.median.mutations,c.intergenic.point.mutations, c.dS.mutations)
ggsave(c.all.mut.plot4,filename='../results/figures/all-mutations4.pdf')

##c.all.mut.plot5 <- plot.aerobic.vs.anaerobic.muts3(c.aerobic.mutations, c.anaerobic.mutations,c.median.mutations,c.intergenic.mutations, c.dS.mutations)
c.all.mut.plot5 <- plot.aerobic.vs.anaerobic.muts3(c.aerobic.mutations, c.anaerobic.mutations,c.median.mutations,c.intergenic.point.mutations, c.dS.mutations)
ggsave(c.all.mut.plot5,filename='../results/figures/all-mutations5.pdf')

#######

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

##########################################################################
## look at accumulation of stars over time for genes in the different proteome
## sectors.
## in other words, look at the rates at which the mutations occur over time.
## plot cumulative sum of anaerobic and aerobic dS and dN in each population.

## get proteome sector assignments from Hui et al. 2015 Supplementary Table 2.
## I saved a reduced version of the data.
proteome.assignments <- read.csv('../data/Hui-2015-proteome-section-assignments.csv',as.is=TRUE)
REL606.proteome.assignments <- inner_join(REL606.genes,proteome.assignments)

## add proteome assignment to mutation.data.
sector.mut.data <- inner_join(mutation.data,REL606.proteome.assignments)

##six sectors:  "A" "S" "O" "U" "R" "C"

A.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='A')
S.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='S')
O.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='O')
U.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='U')
R.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='R')
C.sector.mut.data <- filter(sector.mut.data,Sector.assigned=='C')

A.length <- sum(A.sector.mut.data$length)
S.length <- sum(S.sector.mut.data$length)
O.length <- sum(O.sector.mut.data$length)
U.length <- sum(U.sector.mut.data$length)
R.length <- sum(R.sector.mut.data$length)
C.length <- sum(C.sector.mut.data$length)

c.A.muts <- calc.cumulative.muts(A.sector.mut.data, A.length)
c.S.muts <- calc.cumulative.muts(S.sector.mut.data, S.length)
c.O.muts <- calc.cumulative.muts(O.sector.mut.data, O.length)
c.U.muts <- calc.cumulative.muts(U.sector.mut.data, U.length)
c.R.muts <- calc.cumulative.muts(R.sector.mut.data, R.length)
c.C.muts <- calc.cumulative.muts(C.sector.mut.data, C.length)

plot.sector.muts <- function(A.data,S.data,O.data,U.data,R.data,C.data,all.data,log=TRUE) {
  my.plot <- ggplot(A.data) +
  theme_classic() +
  facet_wrap(.~Population,scales='fixed') +
  ylab('Cumulative number of mutations, normalized by target size') +
  xlab('Generations (x 10,000)')

  if (log) {
    my.plot <- my.plot +
    geom_point(data=A.data,aes(x=Generation,y=log(normalized.cs)),color='black',size=0.2) +
    geom_point(data=S.data,aes(x=Generation,y=log(normalized.cs)),color='red',size=0.2) +
    geom_point(data=O.data,aes(x=Generation,y=log(normalized.cs)),color='blue',size=0.2) +
    geom_point(data=U.data,aes(x=Generation,y=log(normalized.cs)),color='green',size=0.2) +
    geom_point(data=R.data,aes(x=Generation,y=log(normalized.cs)),color='yellow',size=0.2) +
    geom_point(data=all.data,aes(x=Generation,y=log(normalized.cs)),color='grey',size=0.2) +
    geom_point(data=C.data,aes(x=Generation,y=log(normalized.cs)),color='orange',size=0.2)
  } else {
    my.plot <- my.plot +
    geom_point(data=R.data,aes(x=Generation,y=normalized.cs),color='yellow',size=0.2) +
    geom_point(data=C.data,aes(x=Generation,y=normalized.cs),color='orange',size=0.2) +
    geom_point(data=O.data,aes(x=Generation,y=normalized.cs),color='blue',size=0.2) +
    geom_point(data=S.data,aes(x=Generation,y=normalized.cs),color='red',size=0.2) +
    geom_point(data=A.data,aes(x=Generation,y=normalized.cs),color='black',size=0.2) +
    geom_point(data=all.data,aes(x=Generation,y=normalized.cs),color='grey',size=0.2) +
    geom_point(data=U.data,aes(x=Generation,y=normalized.cs),color='green',size=0.2)
  }
}

log.sector.plot <- plot.sector.muts(c.A.muts,c.S.muts,c.O.muts,c.U.muts,c.R.muts,c.C.muts,c.mutations,log=TRUE)
ggsave(log.sector.plot,filename='../results/figures/log-sector-plot.pdf')

sector.plot <- plot.sector.muts(c.A.muts,c.S.muts,c.O.muts,c.U.muts,c.R.muts,c.C.muts,c.mutations,log=FALSE)
ggsave(sector.plot,filename='../results/figures/sector-plot.pdf')

##########################################################################
## look at accumulation of stars over time for genes in different eigengenes
## inferred by Wytock and Motter (2018).
## in other words, look at the rates at which the mutations occur over time.
## plot cumulative sum in each population.

## get eigengene sector assignments from Wytock and Motter (2018) Supplementary File 1.
## I saved a reduced version of the data.

eigengenes <- read.csv('../data/Wytock2018-eigengenes.csv',as.is=TRUE)
REL606.eigengenes <- inner_join(REL606.genes,eigengenes)

## add eigengene assignment to mutation.data.
eigengene.mut.data <- inner_join(mutation.data,REL606.eigengenes)

eigengene1.mut.data <- filter(eigengene.mut.data,Eigengene==1)
eigengene2.mut.data <- filter(eigengene.mut.data,Eigengene==2)
eigengene3.mut.data <- filter(eigengene.mut.data,Eigengene==3)
eigengene4.mut.data <- filter(eigengene.mut.data,Eigengene==4)
eigengene5.mut.data <- filter(eigengene.mut.data,Eigengene==5)
eigengene6.mut.data <- filter(eigengene.mut.data,Eigengene==6)
eigengene7.mut.data <- filter(eigengene.mut.data,Eigengene==7)
eigengene8.mut.data <- filter(eigengene.mut.data,Eigengene==8)
eigengene9.mut.data <- filter(eigengene.mut.data,Eigengene==9)

eigen1.length <- sum(eigengene1.mut.data$length)
eigen2.length <- sum(eigengene2.mut.data$length)
eigen3.length <- sum(eigengene3.mut.data$length)
eigen4.length <- sum(eigengene4.mut.data$length)
eigen5.length <- sum(eigengene5.mut.data$length)
eigen6.length <- sum(eigengene6.mut.data$length)
eigen7.length <- sum(eigengene7.mut.data$length)
eigen8.length <- sum(eigengene8.mut.data$length)
eigen9.length <- sum(eigengene9.mut.data$length)

c.eigen1.muts <- calc.cumulative.muts(eigengene1.mut.data, eigen1.length)
c.eigen2.muts <- calc.cumulative.muts(eigengene2.mut.data, eigen2.length)
c.eigen3.muts <- calc.cumulative.muts(eigengene3.mut.data, eigen3.length)
c.eigen4.muts <- calc.cumulative.muts(eigengene4.mut.data, eigen4.length)
c.eigen5.muts <- calc.cumulative.muts(eigengene5.mut.data, eigen5.length)
c.eigen6.muts <- calc.cumulative.muts(eigengene6.mut.data, eigen6.length)
c.eigen7.muts <- calc.cumulative.muts(eigengene7.mut.data, eigen7.length)
c.eigen8.muts <- calc.cumulative.muts(eigengene8.mut.data, eigen8.length)
c.eigen9.muts <- calc.cumulative.muts(eigengene9.mut.data, eigen9.length)

plot.eigen.muts <- function(eigen1.data,
                            eigen2.data,
                            eigen3.data,
                            eigen4.data,
                            eigen5.data,
                            eigen6.data,
                            eigen7.data,
                            eigen8.data,
                            eigen9.data,
                            all.data,
                            log=TRUE) {
  my.plot <- ggplot(eigen1.data) +
  theme_classic() +
  facet_wrap(.~Population,scales='fixed') +
  ylab('Cumulative number of mutations, normalized by target size') +
  xlab('Generations (x 10,000)')

  if (log) {
    my.plot <- my.plot +
    geom_point(data=eigen1.data,aes(x=Generation,y=log(normalized.cs)),color='red',size=0.2) +
    geom_point(data=eigen2.data,aes(x=Generation,y=log(normalized.cs)),color='orange',size=0.2) +
    geom_point(data=eigen3.data,aes(x=Generation,y=log(normalized.cs)),color='yellow',size=0.2) +
    geom_point(data=eigen4.data,aes(x=Generation,y=log(normalized.cs)),color='green',size=0.2) +
    geom_point(data=eigen5.data,aes(x=Generation,y=log(normalized.cs)),color='cyan',size=0.2) +
    geom_point(data=eigen6.data,aes(x=Generation,y=log(normalized.cs)),color='blue',size=0.2) +
    geom_point(data=eigen7.data,aes(x=Generation,y=log(normalized.cs)),color='violet',size=0.2) +
    geom_point(data=eigen8.data,aes(x=Generation,y=log(normalized.cs)),color='pink',size=0.2) +
    geom_point(data=eigen9.data,aes(x=Generation,y=log(normalized.cs)),color='black',size=0.2) +
    geom_point(data=all.data,aes(x=Generation,y=log(normalized.cs)),color='grey',size=0.2)
  } else {
    my.plot <- my.plot +
    geom_point(data=eigen1.data,aes(x=Generation,y=normalized.cs),color='red',size=0.2) +
    geom_point(data=eigen2.data,aes(x=Generation,y=normalized.cs),color='orange',size=0.2) +
    geom_point(data=eigen3.data,aes(x=Generation,y=normalized.cs),color='yellow',size=0.2) +
    geom_point(data=eigen4.data,aes(x=Generation,y=normalized.cs),color='green',size=0.2) +
    geom_point(data=eigen5.data,aes(x=Generation,y=normalized.cs),color='cyan',size=0.2) +
    geom_point(data=eigen6.data,aes(x=Generation,y=normalized.cs),color='blue',size=0.2) +
    geom_point(data=eigen7.data,aes(x=Generation,y=normalized.cs),color='violet',size=0.2) +
    geom_point(data=eigen8.data,aes(x=Generation,y=normalized.cs),color='pink',size=0.2) +
    geom_point(data=eigen9.data,aes(x=Generation,y=normalized.cs),color='black',size=0.2) +
    geom_point(data=all.data,aes(x=Generation,y=normalized.cs),color='grey',size=0.2)

  }
}

log.eigen.plot <- plot.eigen.muts(c.eigen1.muts,
                                c.eigen2.muts,
                                c.eigen3.muts,
                                c.eigen4.muts,
                                c.eigen5.muts,
                                c.eigen6.muts,
                                c.eigen7.muts,
                                c.eigen8.muts,
                                c.eigen9.muts,
                                c.mutations,log=TRUE)
ggsave(log.eigen.plot,filename='../results/figures/log-eigen-plot.pdf')

eigen.plot <- plot.eigen.muts(c.eigen1.muts,
                                c.eigen2.muts,
                                c.eigen3.muts,
                                c.eigen4.muts,
                                c.eigen5.muts,
                                c.eigen6.muts,
                                c.eigen7.muts,
                                c.eigen8.muts,
                                c.eigen9.muts,
                                c.mutations,log=FALSE)
ggsave(eigen.plot,filename='../results/figures/eigen-plot.pdf')

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
ggsave(log.rando.plot,filename='../results/figures/log-rando-plot.pdf')
rando.plot <- plot.random.subsets(full.mutation.data,log=FALSE)
ggsave(rando.plot,filename='../results/figures/rando-plot.pdf')
