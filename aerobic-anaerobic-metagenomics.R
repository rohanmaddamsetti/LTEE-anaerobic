## aerobic-anaerobic-metagenomics.R by Rohan Maddamsetti.

## IMPORTANT TODO: numbers of aerobic and anaerobic genes don't exactly match those in
## measureTargetSize.py. Debug this before publication.

## The goal of this script is to examine evidence for contingency,
## hypothesized to be caused by ArcAB mutations.

## Hypothesis 1: Nkrumah hypothesizes that anaerobic mutations tend to fix more,
## or that there is increased parallelism in those genes,
## after the fixation of arcA or arcB mutations.

## Hypothesis 2: look at mutations arising before and after arcA fixation.
##        and ask if there is an enrichment afterward of anaerobic-specific
##        genes. (could star CDF get at this answer?)


## Hypothesis 3: first look at identity of genes fixing in the 5000 generation
## window after arcA
## WHEN arcA fixes late (ara+4). THEN ask if the same pattern
## is seen in the same window
## after early. 
   
## 2) Do the same analysis, with arcB, anchoring with Ara+1 (
## as the late occurrence).
## Since the hypothesis is that this mutation is doing the
## same thing as arcA, do the same comparison of parallelism.

## 4) count number of stars over time, and plot.

## 5) examine the mutations in the cohorts.

################################################

library(tidyverse)

## get anaerobic-specific and aerobic-specific genes
## (written out by ArcAnalysisScript.R)
aerobic.genes <- read.csv('../results/aerobic-specific-genes.csv',as.is=TRUE)
anaerobic.genes <- read.csv('../results/anaerobic-specific-genes.csv',as.is=TRUE)

## import genes in the table created by:
## conda activate ltee-metagenomix
## cd LTEE-metagenomic-repo
## python rohan-write-csv.py > ../results/LTEE-metagenome-mutations.csv

mutation.data <- read.csv('../results/LTEE-metagenome-mutations.csv',header=TRUE,as.is=TRUE) %>%
mutate(anaerobic=(Gene %in% anaerobic.genes$gene)) %>%
mutate(aerobic=(Gene %in% aerobic.genes$gene)) %>%
mutate(Generation=t0/10000)

## numbers gotten by running measureTargetSize.py.
## Use these to normalize cumulative mutations over time.
anaerobic.synon.sites <- 103549
anaerobic.nonsynon.sites <- 343951
aerobic.synon.sites <- 50022.66666666666
aerobic.nonsynon.sites <- 168705.33333333337
total.synon.sites <- 890787
total.nonsynon.sites <- 3001647


muts.before.gene.fixation <- function(mutation.data, pop, gene) {
  ## return mutations fixing in the 5000 generation window before 'gene' occurs.
  gene.fixations <- filter(mutation.data,Gene==gene) %>%
  filter(fixation==TRUE)
  my.pop.gene.fixations <-  filter(gene.fixations,Population==pop)
  
  my.pop.gene.origin.times <- my.pop.gene.fixations$t0
  all.mutations.before.gene <- data.frame()
  
  ## include case where there are multiple gene fixations in the pop.
  for (origin.time in my.pop.gene.origin.times) {
    cur.mutations.before.gene <- filter(mutation.data,Population==pop) %>%
    filter(tf>=origin.time-5000) %>%
    filter(tf<=origin.time) %>%
    arrange(tf,t0)
    
    all.mutations.before.gene <- rbind(all.mutations.before.gene,cur.mutations.before.gene)
  }
  return(all.mutations.before.gene)
}

muts.after.gene.fixation <- function(mutation.data, pop, gene) {
  ## return function returning mutations fixing in the 5000 generation window after gene.
  gene.fixations <- filter(mutation.data,Gene==gene) %>%
  filter(fixation==TRUE)
  my.pop.gene.fixations <-  filter(gene.fixations,Population==pop)
  
  my.pop.gene.fixation.times <- my.pop.gene.fixations$tf
  all.mutations.after.gene <- data.frame()
  
  ## include case where there are multiple gene fixations in the pop.
  for (fix.time in my.pop.gene.fixation.times) {
    cur.mutations.after.gene <- filter(mutation.data,Population==pop) %>%
    filter(tf>=fix.time) %>%
    filter(tf<=fix.time+5000) %>%
    arrange(tf,t0)
    
    all.mutations.after.gene <- rbind(all.mutations.after.gene,cur.mutations.after.gene)
  }
  return(all.mutations.after.gene)
}

print.mutations.before.gene <- function(mutation.data, pop, gene) {
  mutations.before.gene <- muts.before.gene.fixation(mutation.data,pop,gene) %>%
  filter(Annotation=='missense')

  anaerobic.loci <- unique(filter(mutations.before.gene,anaerobic==TRUE)$Gene)
  print(paste('anaerobic genes before sweep of gene',gene,'in pop',pop))
  print(length(anaerobic.loci))
  
  aerobic.loci <- unique(filter(mutations.before.gene,aerobic==TRUE)$Gene)
  print(paste('aerobic genes before sweep of gene',gene, 'in pop',pop))
  print(length(aerobic.loci))

  return(mutations.before.gene)
}

print.mutations.after.gene <- function(mutation.data, pop, gene) {
  mutations.after.gene <- muts.after.gene.fixation(mutation.data,pop,gene) %>%
  filter(Annotation=='missense')

  anaerobic.loci <- unique(filter(mutations.after.gene,anaerobic==TRUE)$Gene)
  print(paste('anaerobic genes after sweep of gene',gene,'in pop',pop))
  print(length(anaerobic.loci))
  
  aerobic.loci <- unique(filter(mutations.after.gene,aerobic==TRUE)$Gene)
  print(paste('aerobic genes after sweep of gene',gene,'in pop',pop))
  print(length(aerobic.loci))

  return(mutations.after.gene)
}

####### Look at mutations before/after arcA.
### results for Ara+4
p4.mutations.before.arcA <- print.mutations.before.gene(mutation.data,'Ara+4','arcA')
p4.mutations.after.arcA <- print.mutations.after.gene(mutation.data,'Ara+4','arcA')

## results for Ara+3
p3.mutations.before.arcA <- print.mutations.before.gene(mutation.data,'Ara+3','arcA')
p3.mutations.after.arcA <- print.mutations.after.gene(mutation.data,'Ara+3','arcA')

### results for Ara-6
m6.mutations.before.arcA <- print.mutations.before.gene(mutation.data,'Ara-6','arcA')
m6.mutations.after.arcA <- print.mutations.after.gene(mutation.data,'Ara-6','arcA')

### results for Ara+1
p1.mutations.before.arcA <- print.mutations.before.gene(mutation.data,'Ara+1','arcA')
p1.mutations.after.arcA <- print.mutations.after.gene(mutation.data,'Ara+1','arcA')

### results for Ara+2
p2.mutations.before.arcA <- print.mutations.before.gene(mutation.data,'Ara+1','arcA')
p2.mutations.after.arcA <- print.mutations.after.gene(mutation.data,'Ara+1','arcA')

## bind mutations occurring just BEFORE arcA together in a dataframe.
all.muts.before.arcA <- rbind(m6.mutations.before.arcA,
p1.mutations.before.arcA,
p2.mutations.before.arcA,
p3.mutations.before.arcA,
p4.mutations.before.arcA)

## examine parallelism in these genes across populations.
parallelism.before.arcA <- all.muts.before.arcA %>% group_by(Gene) %>%
summarise(count=n()) %>% filter(count>1) %>% arrange(desc(count))

anaerobic.parallelism.before.arcA <- filter(parallelism.before.arcA,Gene %in% anaerobic.genes$gene)

aerobic.parallelism.before.arcA <- filter(parallelism.before.arcA,Gene %in% aerobic.genes$gene)

anaerobic.parallelism.before.arcA
aerobic.parallelism.before.arcA

## bind mutations occurring just AFTER arcA together in a dataframe.
all.muts.after.arcA <- rbind(m6.mutations.after.arcA,
p1.mutations.after.arcA,
p2.mutations.after.arcA,
p3.mutations.after.arcA,
p4.mutations.after.arcA)

## again, examine parallelism in these genes across populations.
parallelism.after.arcA <- all.muts.after.arcA %>% group_by(Gene) %>%
summarise(count=n()) %>% filter(count>1) %>% arrange(desc(count))

anaerobic.parallelism.after.arcA <- filter(parallelism.after.arcA,Gene %in% anaerobic.genes$gene)

aerobic.parallelism.after.arcA <- filter(parallelism.after.arcA,Gene %in% aerobic.genes$gene)

anaerobic.parallelism.after.arcA
aerobic.parallelism.after.arcA

###########
###########
## now do the same analysis for arcB.

####### Look at mutations before/after arcB.
### results for Ara+1
p1.mutations.before.arcB <- print.mutations.before.gene(mutation.data,'Ara+1','arcB')
p1.mutations.after.arcB <- print.mutations.after.gene(mutation.data,'Ara+1','arcB')

### results for Ara-1
m1.mutations.before.arcB <- print.mutations.before.gene(mutation.data,'Ara-1','arcB')
m1.mutations.after.arcB <- print.mutations.after.gene(mutation.data,'Ara-1','arcB')

### results for Ara-4
m4.mutations.before.arcB <- print.mutations.before.gene(mutation.data,'Ara-4','arcB')
m4.mutations.after.arcB <- print.mutations.after.gene(mutation.data,'Ara-4','arcB')

## results for Ara+3
p3.mutations.before.arcB <- print.mutations.before.gene(mutation.data,'Ara+3','arcB')
p3.mutations.after.arcB <- print.mutations.after.gene(mutation.data,'Ara+3','arcB')

### results for Ara+6
p6.mutations.before.arcB <- print.mutations.before.gene(mutation.data,'Ara+6','arcB')
p6.mutations.after.arcB <- print.mutations.after.gene(mutation.data,'Ara+6','arcB')

##############
all.muts.before.arcB <- rbind(
  p1.mutations.before.arcB,
  m1.mutations.before.arcB,
  m4.mutations.before.arcB,
  p3.mutations.before.arcB,
  p6.mutations.before.arcB)

parallelism.before.arcB <- all.muts.before.arcB %>% group_by(Gene) %>%
summarise(count=n()) %>% filter(count>1) %>% arrange(desc(count))

anaerobic.parallelism.before.arcB <- filter(parallelism.before.arcB,Gene %in% anaerobic.genes$gene)

aerobic.parallelism.before.arcB <- filter(parallelism.before.arcB,Gene %in% aerobic.genes$gene)

anaerobic.parallelism.before.arcB
aerobic.parallelism.before.arcB

##
all.muts.after.arcB <- rbind(
  p1.mutations.after.arcB,
  m1.mutations.after.arcB,
  m4.mutations.after.arcB,
  p3.mutations.after.arcB,
  p6.mutations.after.arcB)


parallelism.after.arcB <- all.muts.after.arcB %>% group_by(Gene) %>%
summarise(count=n()) %>% filter(count>1) %>% arrange(desc(count))

anaerobic.parallelism.after.arcB <- filter(parallelism.after.arcB,Gene %in% anaerobic.genes$gene)

aerobic.parallelism.after.arcB <- filter(parallelism.after.arcB,Gene %in% aerobic.genes$gene) 

anaerobic.parallelism.after.arcB
aerobic.parallelism.after.arcB

## doesn't look like any arcA contingency.
## But perhaps contingency associated with arcB?
## chiA seems like the best candidate.
## looking at the timecourses for this gene,
## seems like some pattern... but hard to say if meaningful.

########################################
########################################
## look at accumulation of stars over time.
## in other words, look at the rates at which the mutations
## occur over time.
## plot cumulative sum of anaerobic and aerobic dS and dN
## in each population.
cumsum.per.pop.helper.func <- function(df) {
  df %>% arrange(t0) %>% group_by(Population,Generation) %>% summarize(count=n()) %>% mutate(cs=cumsum(count))
}

calc.cumulative.muts <- function(data, normalization.constant) {
  data %>% split(.$Population) %>%
  map_dfr(.f=cumsum.per.pop.helper.func) %>%
  mutate(normalized.cs=cs/normalization.constant)
}

dN.mutation.data <- filter(mutation.data,
                           Annotation=='missense')
dS.mutation.data <- filter(mutation.data,
                           Annotation=='synonymous')

extinct.mutation.data <- filter(mutation.data,final_frequency==0)
extinct.dS.mutation.data <- filter(mutation.data,
                                   Annotation=='missense',
                                   final_frequency==0)
extinct.dN.mutation.data <- filter(mutation.data,
                                   Annotation=='synonymous',
                                   final_frequency==0)
extinct.aerobic.mutation.data <- filter(extinct.mutation.data,aerobic==TRUE)
extinct.anaerobic.mutation.data <- filter(extinct.mutation.data,anaerobic==TRUE)
extinct.aerobic.dN.mutation.data <- filter(extinct.dN.mutation.data,aerobic==TRUE)
extinct.anaerobic.dN.mutation.data <- filter(extinct.dN.mutation.data,anaerobic==TRUE)
extinct.aerobic.dS.mutation.data <- filter(extinct.dS.mutation.data,aerobic==TRUE)
extinct.anaerobic.dS.mutation.data <- filter(extinct.dS.mutation.data,anaerobic==TRUE)

c.extinct.aerobic.mutations <- calc.cumulative.muts(extinct.aerobic.mutation.data,
                                                             aerobic.gene.length)
c.extinct.anaerobic.mutations <- calc.cumulative.muts(extinct.anaerobic.mutation.data,
                                                             anaerobic.gene.length)
c.extinct.aerobic.dN.mutations <- calc.cumulative.muts(extinct.aerobic.dN.mutation.data,
                                                             total.nonsynon.sites)
c.extinct.anaerobic.dN.mutations <- calc.cumulative.muts(extinct.anaerobic.dN.mutation.data,
                                                             total.nonsynon.sites)
c.extinct.aerobic.dS.mutations <- calc.cumulative.muts(extinct.aerobic.dS.mutation.data,
                                                             total.synon.sites)
c.extinct.anaerobic.dS.mutations <- calc.cumulative.muts(extinct.anaerobic.dS.mutation.data,
                                                             total.synon.sites)

c.dN.mutations <- calc.cumulative.muts(dN.mutation.data,total.nonsynon.sites)

low.freq.dS.mutation.data <- filter(dS.mutation.data,final_frequency<0.2)
med.freq.dS.mutation.data <- filter(dS.mutation.data,final_frequency>0.2) %>%
filter(final_frequency<0.8)
high.freq.dS.mutation.data <- filter(dS.mutation.data,final_frequency>0.8)

c.dN.mutations <- calc.cumulative.muts(dN.mutation.data,total.nonsynon.sites)
c.dS.mutations <- calc.cumulative.muts(dS.mutation.data,total.synon.sites)

c.low.freq.dS.mutations <- calc.cumulative.muts(low.freq.dS.mutation.data,
                                                total.synon.sites)
c.med.freq.dS.mutations <- calc.cumulative.muts(med.freq.dS.mutation.data,
                                                total.synon.sites)
c.high.freq.dS.mutations <- calc.cumulative.muts(high.freq.dS.mutation.data,
                                                 total.synon.sites) 

aerobic.low.freq.dS.mutation.data <- filter(low.freq.dS.mutation.data,aerobic==TRUE)
aerobic.med.freq.dS.mutation.data <- filter(med.freq.dS.mutation.data,aerobic==TRUE)
aerobic.high.freq.dS.mutation.data <- filter(high.freq.dS.mutation.data,aerobic==TRUE)

c.aerobic.low.freq.dS.mutations <- calc.cumulative.muts(aerobic.low.freq.dS.mutation.data,
                                                        total.synon.sites)
c.aerobic.med.freq.dS.mutations <- calc.cumulative.muts(aerobic.med.freq.dS.mutation.data,
                                                        total.synon.sites)
c.aerobic.high.freq.dS.mutations <- calc.cumulative.muts(aerobic.high.freq.dS.mutation.data,
                                                         total.synon.sites) 

anaerobic.low.freq.dS.mutation.data <- filter(low.freq.dS.mutation.data,
                                              anaerobic==TRUE)
anaerobic.med.freq.dS.mutation.data <- filter(med.freq.dS.mutation.data,
                                              anaerobic==TRUE)
anaerobic.high.freq.dS.mutation.data <- filter(high.freq.dS.mutation.data,
                                               anaerobic==TRUE)

c.anaerobic.low.freq.dS.mutations <- calc.cumulative.muts(anaerobic.low.freq.dS.mutation.data,
                                                          total.synon.sites)
c.anaerobic.med.freq.dS.mutations <- calc.cumulative.muts(anaerobic.med.freq.dS.mutation.data,
                                                          total.synon.sites)
c.anaerobic.high.freq.dS.mutations <- calc.cumulative.muts(anaerobic.high.freq.dS.mutation.data,
                                                           total.synon.sites) 


aerobic.mutation.data <- filter(mutation.data,aerobic==TRUE)
aerobic.dN.mutation.data <- filter(dN.mutation.data,aerobic==TRUE)
aerobic.dS.mutation.data <- filter(dS.mutation.data,aerobic==TRUE)

c.aerobic.dN.mutations <- calc.cumulative.muts(aerobic.dN.mutation.data,
                                               aerobic.nonsynon.sites)

c.aerobic.dS.mutations <- calc.cumulative.muts(aerobic.dS.mutation.data,
                                               aerobic.synon.sites)

anaerobic.mutation.data <- filter(mutation.data,anaerobic==TRUE)
anaerobic.dN.mutation.data <- filter(dN.mutation.data,anaerobic==TRUE)
anaerobic.dS.mutation.data <- filter(dS.mutation.data,anaerobic==TRUE)

c.anaerobic.dN.mutations <- calc.cumulative.muts(anaerobic.dN.mutation.data,
                                                 anaerobic.nonsynon.sites)

c.anaerobic.dS.mutations <- calc.cumulative.muts(anaerobic.dS.mutation.data,
                                                 anaerobic.synon.sites)

## plot cumulative plot of anaerobic dN and aerobic dN on top of each other.
c.dN.plot <- ggplot(c.aerobic.dN.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.anaerobic.dN.mutations,aes(x=Generation,y=normalized.cs,color='red'))
ggsave(c.dN.plot,filename='../results/figures/just-dN.pdf')

c.dS.plot <- ggplot(c.aerobic.dS.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.anaerobic.dS.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.dS.plot,filename='../results/figures/just-dS.pdf')

## This is weird! Check out if pattern holds if I split up by final frequency at 60K generations.
## The curves do seem to change depending on final frequency.

c.low.freq.dS.plot <- ggplot(c.aerobic.low.freq.dS.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.anaerobic.low.freq.dS.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.low.freq.dS.plot,filename='../results/figures/lowfreq-dS.pdf')

c.med.freq.dS.plot <- ggplot(c.aerobic.med.freq.dS.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.anaerobic.med.freq.dS.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.med.freq.dS.plot,filename='../results/figures/medfreq-dS.pdf')

c.high.freq.dS.plot <- ggplot(c.aerobic.high.freq.dS.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.anaerobic.high.freq.dS.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.high.freq.dS.plot,filename='../results/figures/highfreq-dS.pdf')




c.total.dN.dS.plot <- ggplot(c.dN.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.dS.mutations,aes(x=Generation,y=normalized.cs,color='dS'))

## based on redoing the binomial test on genomes, it seems the difference in result
## is NOT driven by target size. rather, it is whether or not all kinds of mutations
## are included, and not just restricting to point mutations.

## from ArcAnalysisScript.R:
## anaerobic-specific gene length: 457593
## aerobic-specific gene length: 219286
anaerobic.gene.length <- 457593
aerobic.gene.length <- 219286

c.anaerobic.mutations.no.dS <- calc.cumulative.muts(anaerobic.mutation.data,
                                                    anaerobic.gene.length)

c.aerobic.mutations.no.dS <- calc.cumulative.muts(aerobic.mutation.data,
                                                  aerobic.gene.length)

c.no.dS.plot <- ggplot(c.aerobic.mutations.no.dS,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
xlabel
geom_point(data=c.anaerobic.mutations.no.dS,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.no.dS.plot,filename='../results/figures/all-mutations-but-dS.pdf')

c.aerobic.mutations <- calc.cumulative.muts(aerobic.mutation.data,
                                            aerobic.gene.length)
c.anaerobic.mutations <- calc.cumulative.muts(anaerobic.mutation.data,
                                              anaerobic.gene.length)
c.all.mut.plot <- ggplot(c.aerobic.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.anaerobic.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.all.mut.plot,filename='../results/figures/all-mutations.pdf')

## does this plot depend on final frequency?
high.freq.aerobic.muts <- filter(aerobic.mutation.data,final_frequency>0.8)
high.freq.anaerobic.muts <- filter(anaerobic.mutation.data,final_frequency>0.8)

c.high.freq.aerobic.mutations <- calc.cumulative.muts(high.freq.aerobic.muts,
                                            aerobic.gene.length)
c.high.freq.anaerobic.mutations <- calc.cumulative.muts(high.freq.anaerobic.muts,
                                              anaerobic.gene.length)

c.high.freq.all.mut.plot <- ggplot(c.high.freq.aerobic.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.high.freq.anaerobic.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.high.freq.all.mut.plot,filename='../results/figures/highfreq-all-muts.pdf')

low.freq.aerobic.muts <- filter(aerobic.mutation.data,final_frequency<0.2)
low.freq.anaerobic.muts <- filter(anaerobic.mutation.data,final_frequency<0.2)

c.low.freq.aerobic.mutations <- calc.cumulative.muts(low.freq.aerobic.muts,
                                            aerobic.gene.length)
c.low.freq.anaerobic.mutations <- calc.cumulative.muts(low.freq.anaerobic.muts,
                                              anaerobic.gene.length)

c.low.freq.all.mut.plot <- ggplot(c.low.freq.aerobic.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.low.freq.anaerobic.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.low.freq.all.mut.plot,filename='../results/figures/lowfreq-all-muts.pdf')


### let's plot just extinct mutations. Any insights from this?
c.extinct.mut.plot <- ggplot(c.extinct.aerobic.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.extinct.anaerobic.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.extinct.mut.plot,filename='../results/figures/extinct-all-muts.pdf')

c.extinct.dN.plot <- ggplot(c.extinct.aerobic.dN.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.extinct.anaerobic.dN.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.extinct.dN.plot,filename='../results/figures/extinct-dN.pdf')

c.extinct.dS.plot <- ggplot(c.extinct.aerobic.dS.mutations,aes(x=Generation,y=normalized.cs)) +
theme_classic() +
geom_point() +
facet_wrap(.~Population,scales='free') +
geom_point(data=c.extinct.anaerobic.dS.mutations,aes(x=Generation,y=normalized.cs,color='anaerobic'))
ggsave(c.extinct.dN.plot,filename='../results/figures/extinct-dS.pdf')




################################################################
## Doesn't look like anaerobic or aerobic mutations cluster in the same sweep.

## let's examine cohorts of mutations as a check.
## try sorting mutations into cohorts based on Population, tf, and t0.
## could consider clustering and visualizing cohorts or something.

## IMPORTANT NOTE: these are really just mutations occurring at the same time,
## and not necessarily true 'cohorts' on the same lineage
## (though those will count here too).
## They may not actually be in the same background.
## For instance, clonally interferening arcB mutations would count here.

sorted.mutation.data <- mutation.data %>% #group_by(Population) %>%
select(-Allele) %>% arrange(Population,tf, t0)

print.cohorts <- function(sorted.mutation.data) {
  for (p in unique(sorted.mutation.data$Population)) {
    pop.mutations <- filter(sorted.mutation.data,Population==p)
    for (final.time in unique(pop.mutations$tf)) {
      final.time.mutations <- filter(pop.mutations,tf==final.time)
      for (initial.time in unique(final.time.mutations$t0)) {
        cohort <- filter(final.time.mutations,t0==initial.time)
        print(paste('t0:',unique(cohort$t0),'tf:',unique(cohort$tf),'cohort:'))
        print(cohort)
      }
    }
  } 
}

## summarize cohorts, and sort by fraction of anaerobic dN.
cohort.summary <- sorted.mutation.data %>%
filter(Annotation=='missense') %>%
group_by(Population, tf, t0) %>%
summarize(anaerobic.count=sum(anaerobic),
          aerobic.count=sum(aerobic),
          total.count=n()) %>%
filter(tf<1000000) %>% ## ignore coexistence.
mutate(anaerobic.frac=anaerobic.count/total.count) %>%
filter(total.count>1) %>%
arrange(desc(anaerobic.frac),desc(total.count))

anaerobic.cohorts <- filter(cohort.summary,anaerobic.count>=1)

## let's do a quick binomial test to see if anaerobic genes
## co-occur more often than expected.
length(aerobic.genes$gene) ## number of anerobic-specific genes: 227
length(anaerobic.genes$gene) ## number of anerobic-specific genes: 345
## total number of mutated genes: 3970
length(unique(sorted.mutation.data$Gene))

345/3970 ## 0.0869 ## rough P(hit anaerobic gene by random mutation)

## 2,291 cohorts examined. multiply for a bonferroni correction.
2291*dbinom(2, size=2, prob=(345/3970)) ## 17.3 hits expected
2291*pbinom(2, size=3, prob=(345/3970),lower.tail=F) ## 1.5 hits expected
2291*pbinom(3, size=5, prob=(345/3970),lower.tail=F) ## 0.6 hits expected
2291*pbinom(4, size=7, prob=(345/3970),lower.tail=F) ## 0.2 hits expected
2291*pbinom(6, size=12, prob=(345/3970),lower.tail=F) ## 0.0457 hits expected
2291*pbinom(5, size=10, prob=(345/3970),lower.tail=F) ## 0.15 hits expected


## Is this a signal of co-occurrence of anaerobic dN more often that expected?
## This signal only occurs in mutators though,
## when I compare to co-occurring mutations in non.mutators, I don't see anything:

non.mutators <- c('Ara-5','Ara+2','Ara+4','Ara+1','Ara+5','Ara-6')
non.mutator.cohorts <- filter(cohort.summary,Population %in% non.mutators)
non.mutator.cohorts

## calculate the probability of a given configuration of parallel mutations
## in the metagenomes across populations.
## This seems to be a neat test for positive selection!
## Seems to give the same answer as results in Tenaillon Nature paper,
## and unfortunately not effective for historical contingency per se.

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
