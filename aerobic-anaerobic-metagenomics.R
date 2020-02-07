## aerobic-anaerobic-metagenomics.R by Rohan Maddamsetti.

## IMPORTANT TODO: numbers of aerobic and anaerobic genes don't exactly match those in
## measureTargetSize.py. Debug this before publication.

## IMPORTANT: run anaerobic-anaerobic-genomics.R first, to generate
## anaerobic-specific-genes.csv and aerobic-specific-genes.csv.

## Basic premise.
## count the cumulative number of stars over time, and plot.
## examine different kinds of mutations and genes.

## MAKE SURE THAT ZEROS IN MUTATION COUNTS ARE NOT THROWN OUT BY DPLYR!

## use harmonic average of p-value for aggregation? Or Fisher's method?
## See Daniel Wilson's recent papers.

library(tidyverse)
##########################################################################

## Order nonmutator pops, then hypermutator pops by converting Population to
## factor type and setting the levels.
nonmutator.pops <- c("Ara-5", "Ara-6", "Ara+1", "Ara+2", "Ara+4", "Ara+5")
hypermutator.pops <- c("Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara+3", "Ara+6")


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
    mutate(anaerobic=(Gene %in% anaerobic.genes$Gene)) %>%
    mutate(aerobic=(Gene %in% aerobic.genes$Gene)) %>%
    ## Order nonmutator pops, then hypermutator pops by converting Population to
    ## factor type and setting the levels.
    mutate(Population=factor(Population,levels=c(nonmutator.pops,hypermutator.pops)))

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
## Needed to normalize cumulative mutations in non-coding regions.
##intergenic.length <- 487863

## from aerobic-anaerobic-genomics.R: ########
######################NOTE: DOUBLE CHECK CONSISTENT WITH measureTargetSize.py. output!!!!!

##anaerobic.gene.length <- 457593
##aerobic.gene.length <- 219286
##########################################################################

## Normalization constant calculations.
## IMPORTANT-- THIS STEP IS REALLY IMPORTANT.
## NOT CLEAR HOW TO NORMALIZE IN ORDER TO COMPARE DIFFERENT CLASSES OF MUTATIONS
## "APPLES TO APPLES".


## IMPORTANT TODO: check difference between normalizing by gene length and normalizing
## by synonymous/nonsynonymous opportunities. The end result should be very similar.
## Then consider refactoring to just use REL606_IDs.csv to get gene lengths.

## IMPORTANT TODO: split nonsynonymous into opportunities for missense and nonsense
## mutations, as separate classes, to fit annotation given in Ben Good's data.

## Example of how to try this normalization.
##dN.normalization.const <- total.length * total.nonsynon.sites/(total.synon.sites+total.nonsynon.sites)
##dS.normalization.const <- total.length * total.synon.sites/(total.synon.sites+total.nonsynon.sites)

#####################################################
## Think about cutting this stuff-- internal consistency is most important in the
## comparison to the bootstrapped null.

## numbers gotten by running measureTargetSize.py.
## Use these to normalize cumulative mutations over time.
##target.size.numbers <- read.csv('../results/target_size.csv',header=TRUE,as.is=TRUE)

##anaerobic.length <- filter(target.size.numbers,set=='anaerobic')$total_gene_length
##anaerobic.synon.sites <- filter(target.size.numbers,set=='anaerobic')$synon_sites
##anaerobic.nonsynon.sites <- filter(target.size.numbers,set=='anaerobic')$non_synon_sites

##aerobic.length <- filter(target.size.numbers,set=='aerobic')$total_gene_length
##aerobic.synon.sites <- filter(target.size.numbers,set=='aerobic')$synon_sites
##aerobic.nonsynon.sites <- filter(target.size.numbers,set=='aerobic')$non_synon_sites


##total.length <- filter(target.size.numbers,set=='genome')$total_gene_length
##total.synon.sites <- filter(target.size.numbers,set=='genome')$synon_sites
##total.nonsynon.sites <- filter(target.size.numbers,set=='genome')$non_synon_sites

########################################
## look at accumulation of stars over time.
## in other words, look at the rates at which the mutations occur over time.
## plot cumulative sum of anaerobic and aerobic dS and dN in each population.

cumsum.per.pop.helper.func <- function(df) {
    df %>%
        arrange(t0) %>%
        ## very important: don't drop empty groups because we want to keep zeros.
        group_by(Population,Generation,.drop=FALSE) %>%
        summarize(count=n()) %>%
        mutate(cs=cumsum(count)) %>%
        ungroup()
}

calc.cumulative.muts <- function(data, normalization.constant=NA) {

    ## if normalization.constant is not provided, then
    ## calculate based on gene length by default.
    if (is.na(normalization.constant)) {
        my.genes <- data %>% select(Gene,gene_length) %>% distinct()
        normalization.constant <- sum(my.genes$gene_length)
    }
    
    c.dat <- data %>%
        split(.$Population) %>%
        map_dfr(.f=cumsum.per.pop.helper.func) %>%
        mutate(normalized.cs=cs/normalization.constant) %>%
        ## remove any NA values.
        na.omit()
    return(c.dat)
}
#########################################################################
## calculate cumulative numbers of mutations in each category.
## Then, make plots to see if any interesting patterns emerge.

## for vanilla plotting.
plot.cumulative.muts <- function(mut.data,logscale=TRUE, my.color="black") {
    if (logscale) {
        p <- ggplot(mut.data,aes(x=Generation,y=log10(normalized.cs))) +
            ##ylim(-7,-2) +
            ylab('log[Cumulative number of mutations (normalized)]')
    } else {
        p <- ggplot(mut.data,aes(x=Generation,y=normalized.cs)) +
            ##ylim(0,0.003) +
            ylab('Cumulative number of mutations (normalized)')
    }
    p <- p +
        theme_classic() +
        geom_point(size=0.2, color=my.color) +
        geom_step(size=0.2, color=my.color) +
        facet_wrap(.~Population,scales='free',nrow=4) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3) +
        theme(axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14))
    return(p)
}      

## This plot visualizes a two-tailed test (alpha = 0.05)
## against a bootstrapped null distribution.
## This is what we want to use for publication.
plot.base.layer <- function(data, subset.size=300, N=1000, alpha = 0.05, logscale=FALSE, normalization.constant=NA) {

    ## This function takes the index for the current draw, and samples the data,
    ## generating a random gene set for which to calculate cumulative mutations.
    generate.cumulative.mut.subset <- function(idx) {
        rando.genes <- sample(unique(data$Gene),subset.size)
        mut.subset <- filter(data,Gene %in% rando.genes)
        c.mut.subset <- calc.cumulative.muts(mut.subset, normalization.constant) %>%
            mutate(bootstrap_replicate=idx)
        return(c.mut.subset)
    }

    ## make a dataframe of bootstrapped trajectories.
    ## look at accumulation of stars over time for random subsets of genes.
    ## I want to plot the distribution of cumulative mutations over time for
    ## say, 1000 or 10000 random subsets of genes.

    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),.f=generate.cumulative.mut.subset)

    ## filter out the top alpha/2 and bottom alpha/2 trajectories from each population,
    ## for a two-sided test. default is alpha == 0.05.
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() 
    
    top.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        top_frac(alpha/2) %>%
        select(-final.norm.cs) %>%
        mutate(in.top=TRUE)
    
    bottom.trajectories <- trajectory.summary %>%
        group_by(Population) %>%
        top_frac(-alpha/2) %>%
        select(-final.norm.cs) %>%
        mutate(in.bottom=TRUE)
    
    filtered.trajectories <- bootstrapped.trajectories %>%
        left_join(top.trajectories) %>%
        left_join(bottom.trajectories) %>%
        filter(is.na(in.top)) %>%
        filter(is.na(in.bottom)) %>%
        select(-in.top,-in.bottom)

    if (logscale) {
        p <- ggplot(filtered.trajectories,aes(x=Generation,y=log10(normalized.cs))) +
            ylab('log[Cumulative number of mutations (normalized)]')
    } else {
        p <- ggplot(filtered.trajectories,aes(x=Generation,y=normalized.cs)) +
            ylab('Cumulative number of mutations (normalized)')
    }
    p <- p +
        theme_classic() +
        geom_point(size=0.2, color='gray') +
        facet_wrap(.~Population,scales='free',nrow=4) +
        xlab('Generations (x 10,000)') +
        xlim(0,6.3) +
        theme(axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.x  = element_text(size=14),
              axis.text.y  = element_text(size=14))
    return(p)                
}

## take a ggplot object output by plot.cumulative.muts, and add an extra layer.
add.cumulative.mut.layer <- function(p, layer.df, my.color, logscale=FALSE) {
    if (logscale) {
        p <- p +
            geom_point(data=layer.df, aes(x=Generation,y=log10(normalized.cs)), color=my.color, size=0.2) +
            geom_step(data=layer.df, aes(x=Generation,y=log10(normalized.cs)), size=0.2, color=my.color)
        } else {
            p <- p +
                geom_point(data=layer.df, aes(x=Generation,y=normalized.cs), color=my.color, size=0.2) +
                geom_step(data=layer.df, aes(x=Generation,y=normalized.cs), size=0.2, color=my.color)
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

## let's look at noncoding mutations.
gene.noncoding.mutation.data <- filter(gene.mutation.data,Annotation=='noncoding')


## normalizing constants need to be consistent with the null distributions!
## gene length normalization now calculated by default in calc.cumulative.muts.

## let's look at all mutations within genes.
c.mutations <- calc.cumulative.muts(gene.mutation.data)
                                    
aerobic.mutation.data <- filter(gene.mutation.data,aerobic==TRUE)
anaerobic.mutation.data <- filter(gene.mutation.data,anaerobic==TRUE)

aerobic.dN.mutation.data <- filter(aerobic.mutation.data,
                                   Annotation=='missense')
anaerobic.dN.mutation.data <- filter(anaerobic.mutation.data,
                                     Annotation=='missense')

aerobic.dS.mutation.data <- filter(aerobic.mutation.data,
                                   Annotation=='synonymous')
anaerobic.dS.mutation.data <- filter(anaerobic.mutation.data,
                                     Annotation=='synonymous')

aerobic.nonsense.mutation.data <- filter(aerobic.mutation.data,
                                   Annotation=='nonsense')
anaerobic.nonsense.mutation.data <- filter(anaerobic.mutation.data,
                                     Annotation=='nonsense')

## NOTE APPARENT PARALLELISM IN NON-CODING MUTATIONS!
## INVESTIGATE THIS FURTHER AT THE END.
aerobic.noncoding.mutation.data <- filter(aerobic.mutation.data,
                                      Annotation=='noncoding')
anaerobic.noncoding.mutation.data <- filter(anaerobic.mutation.data,
                                      Annotation=='noncoding')

c.aerobic.mutations <- calc.cumulative.muts(aerobic.mutation.data)
c.anaerobic.mutations <- calc.cumulative.muts(anaerobic.mutation.data)

c.aerobic.dN.mutations <- calc.cumulative.muts(aerobic.dN.mutation.data)
c.anaerobic.dN.mutations <- calc.cumulative.muts(anaerobic.dN.mutation.data)

c.aerobic.dS.mutations <- calc.cumulative.muts(aerobic.dS.mutation.data)
c.anaerobic.dS.mutations <- calc.cumulative.muts(anaerobic.dS.mutation.data)

c.aerobic.nonsense.mutations <- calc.cumulative.muts(aerobic.nonsense.mutation.data)
c.anaerobic.nonsense.mutations <- calc.cumulative.muts(anaerobic.nonsense.mutation.data)

## normalize noncoding mutations by 1 across the board.
c.aerobic.noncoding.mutations <- calc.cumulative.muts(aerobic.noncoding.mutation.data,
                                                      normalization.constant=1)

c.anaerobic.noncoding.mutations <- calc.cumulative.muts(anaerobic.noncoding.mutation.data,
                                                        normalization.constant=1)

#############################################################################
## Figures. Plot real data on top of the random expectations to plot hypothesis tests.

## plot all classes of mutations in aerobic versus anaerobic genes.
## Figure for All Mutation Classes.
c.aerobic.vs.anaerobic.plot <- plot.base.layer(gene.mutation.data) %>%
    add.cumulative.mut.layer(c.aerobic.mutations,my.color='black') %>%
    add.cumulative.mut.layer(c.anaerobic.mutations, my.color='red')
ggsave(filename="../results/figures/AllMutFig2.pdf", plot=c.aerobic.vs.anaerobic.plot)

## make the same plot, for structural variation, indels, nonsense mutations.
sv.indel.nonsense.muts <- gene.mutation.data %>%
    filter(Annotation %in% c("sv", "indel", "nonsense"))
c.sv.indel.nonsense.muts <- calc.cumulative.muts(sv.indel.nonsense.muts)

## aerobic in black
aerobic.sv.indel.nonsense.muts <- sv.indel.nonsense.muts %>%
    filter(aerobic==TRUE)
c.aerobic.sv.indel.nonsense.muts <- calc.cumulative.muts(aerobic.sv.indel.nonsense.muts)

## anaerobic in red
anaerobic.sv.indel.nonsense.muts <- sv.indel.nonsense.muts %>%
    filter(anaerobic==TRUE)
c.anaerobic.sv.indel.nonsense.muts <- calc.cumulative.muts(anaerobic.sv.indel.nonsense.muts)


aerobic.anaerobic.purifying.plot <- plot.base.layer(sv.indel.nonsense.muts) %>%
    add.cumulative.mut.layer(c.aerobic.sv.indel.nonsense.muts, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.sv.indel.nonsense.muts, my.color="red")
ggsave(filename="../results/figures/SVIndelNonsense.pdf", plot=aerobic.anaerobic.purifying.plot)

## plot just dN.
aerobic.anaerobic.dN.plot <- plot.base.layer(gene.dN.mutation.data) %>%
        add.cumulative.mut.layer(c.aerobic.dN.mutations, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.dN.mutations, my.color="red")    
ggsave(filename="../results/figures/dN.pdf", plot=aerobic.anaerobic.dN.plot)


## plot just dS.
aerobic.anaerobic.dS.plot <- plot.base.layer(gene.dS.mutation.data) %>%
    add.cumulative.mut.layer(c.aerobic.dS.mutations, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.dS.mutations, my.color="red") 
ggsave(filename="../results/figures/dS.pdf", plot=aerobic.anaerobic.dS.plot)

## plot just nonsense mutations.
aerobic.anaerobic.nonsense.plot <- plot.base.layer(gene.nonsense.mutation.data) %>%
    add.cumulative.mut.layer(c.aerobic.nonsense.mutations, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.nonsense.mutations, my.color="red") 
ggsave(filename="../results/figures/nonsense.pdf", plot=aerobic.anaerobic.nonsense.plot)

## plot just non-coding mutations.
## normalize by 1 since dividing by gene length doesn't make sense in this case.
## BUG: NORMALIZATION IS STILL OFF!
aerobic.anaerobic.noncoding.plot <- plot.base.layer(gene.noncoding.mutation.data,
                                                    normalization.constant=1) %>%
    add.cumulative.mut.layer(c.aerobic.noncoding.mutations, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.noncoding.mutations, my.color="red")
aerobic.anaerobic.noncoding.plot
ggsave(filename="../results/figures/noncoding.pdf", plot=aerobic.anaerobic.noncoding.plot)

########################################################################################
## LOOK AT PARALLEL EVOLUTION IN NON-CODING MUTATIONS!
## LOOK AT mntH and hupA
aerobic.noncoding.summary <- aerobic.noncoding.mutation.data %>% group_by(Gene,Position) %>%
    summarize(count=n()) %>% arrange(desc(count,Position))

aerobic.noncoding.summary2 <- aerobic.noncoding.mutation.data %>% group_by(Gene) %>%
    summarize(count=n()) %>% filter(count>1) %>% arrange(desc(count,Position))

aerobic.noncoding.summary3 <- filter(aerobic.noncoding.summary, Gene %in% aerobic.noncoding.summary2$Gene)

anaerobic.noncoding.summary <- anaerobic.noncoding.mutation.data %>% group_by(Gene,Position) %>% summarize(count=n()) %>% arrange(desc(count,Position))

anaerobic.noncoding.summary2 <- anaerobic.noncoding.mutation.data %>% group_by(Gene) %>%
    summarize(count=n()) %>% filter(count > 1) %>% arrange(desc(count,Position))

anaerobic.noncoding.summary3 <- filter(anaerobic.noncoding.summary, Gene %in% anaerobic.noncoding.summary2$Gene)

## let's look at all noncoding mutations.
noncoding.mutation.data <- mutation.data %>% filter(Annotation=='noncoding')

noncoding.summary1 <- noncoding.mutation.data %>% group_by(Gene,Position) %>%
    summarize(count=n()) %>% arrange(desc(count,Position))

noncoding.summary2 <- noncoding.mutation.data %>% group_by(Gene) %>%
    summarize(count=n()) %>% filter(count > 3) %>% arrange(desc(count,Position))

noncoding.summary3 <- filter(noncoding.summary1, Gene %in% noncoding.summary2$Gene) 

noncoding.summary4 <- filter(noncoding.summary3,Gene != "intergenic")
noncoding.summary4

mntH.muts <- filter(mutation.data,Gene == 'mntH')
