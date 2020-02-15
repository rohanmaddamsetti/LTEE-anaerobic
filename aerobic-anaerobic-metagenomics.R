## aerobic-anaerobic-metagenomics.R by Rohan Maddamsetti.

## IMPORTANT: numbers of aerobic and anaerobic genes don't exactly match those in
## measureTargetSize.py. At present that script is not used in this analysis at all.

## IMPORTANT: run anaerobic-anaerobic-genomics.R first, to generate
## anaerobic-specific-genes.csv and aerobic-specific-genes.csv.

## Basic premise.
## count the cumulative number of stars over time, and plot.
## examine different kinds of mutations and genes.

## MAKE SURE THAT ZEROS IN MUTATION COUNTS ARE NOT THROWN OUT BY DPLYR!

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

########################################
## function for plotting better y-axis labels.
## see solution here for nice scientific notation on axes.
## https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
fancy_scientific <- function(x) {
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

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

## calculate the tail probabilities of the true cumulative mutation trajectory
## of a given vector of genes (a 'module'), based on resampling
## random sets of genes. Returns both upper tail of null distribution,
## or P(random trajectory >= the actual trajectory).
## Output: a dataframe with three columns: Population, count, p.val
calculate.trajectory.tail.probs <- function(data, gene.vec, N=10000, normalization.constant=NA) {

    ## resamples have the same cardinality as the gene.vec.
    subset.size <- length(gene.vec)
    
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

    bootstrapped.trajectories <- map_dfr(.x=seq_len(N),.f=generate.cumulative.mut.subset)

    gene.vec.data <- data %>% filter(Gene %in% gene.vec)
    data.trajectory <- calc.cumulative.muts(gene.vec.data,normalization.constant)
    data.trajectory.summary <- data.trajectory %>%
        group_by(Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup()
    
    trajectory.summary <- bootstrapped.trajectories %>%
        ## important: don't drop empty groups.
        group_by(bootstrap_replicate, Population,.drop=FALSE) %>%
        summarize(final.norm.cs=max(normalized.cs)) %>%
        ungroup() 

    trajectory.filter.helper <- function(pop.trajectories) {
        pop <- unique(pop.trajectories$Population)
        data.traj <- filter(data.trajectory.summary,Population == pop)
        final.data.norm.cs <- unique(data.traj$final.norm.cs)
        tail.trajectories <- filter(pop.trajectories, final.norm.cs >= final.data.norm.cs)
        return(tail.trajectories)
    }
    
    ## split by Population, then filter for bootstraps > data trajectory.
    uppertail.probs <- trajectory.summary %>%
        split(.$Population) %>%
        map_dfr(.f=trajectory.filter.helper) %>%
        group_by(Population) %>%
        summarize(count=n()) %>%
        mutate(p.val=count/N)
        
    return(uppertail.probs)
}

#########################################################################
## calculate cumulative numbers of mutations in each category.
## Then, make plots to see if any interesting patterns emerge.

## for vanilla plotting.
plot.cumulative.muts <- function(mut.data,logscale=TRUE, my.color="black") {
    if (logscale) {
        p <- ggplot(mut.data,aes(x=Generation,y=log10(normalized.cs))) +
            ylab('log[Cumulative number of mutations (normalized)]')
    } else {
        p <- ggplot(mut.data,aes(x=Generation,y=normalized.cs)) +
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
              axis.text.y  = element_text(size=14)) +
        scale_y_continuous(labels=fancy_scientific)
    
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

#############################################################################
## Analysis and figures.
## Plot real data on top of the random expectations to plot hypothesis tests.

## look at all mutations within genes.
## calculate p-values for all mutations.
aerobic.allmut.pvals <- calculate.trajectory.tail.probs(gene.mutation.data, aerobic.genes$Gene)
anaerobic.allmut.pvals <- calculate.trajectory.tail.probs(gene.mutation.data, anaerobic.genes$Gene)

## plot all classes of mutations in aerobic versus anaerobic genes.
## profile code to see how number of bootstraps affects runtime.
start.time <- Sys.time()

c.aerobic.mutations <- gene.mutation.data %>%
    filter(aerobic==TRUE) %>%
    calc.cumulative.muts()

c.anaerobic.mutations <- gene.mutation.data %>%
    filter(anaerobic==TRUE) %>%
    calc.cumulative.muts()

c.aerobic.vs.anaerobic.plot <- plot.base.layer(gene.mutation.data) %>%
    add.cumulative.mut.layer(c.aerobic.mutations,my.color='black') %>%
    add.cumulative.mut.layer(c.anaerobic.mutations, my.color='red')
ggsave(filename="../results/figures/AllMutFig2.pdf", plot=c.aerobic.vs.anaerobic.plot)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
## end of profiling.
## 7 min for 10,000 bootstraps.
## PDF image is 28MB on disk.

######
## indels, structural variation (mobile elements), and nonsense mutations.

sv.indel.nonsense.muts <- gene.mutation.data %>%
    filter(Annotation %in% c("sv", "indel", "nonsense"))

## calculate p-values.
aerobic.sv.indel.nonsense.pvals <- calculate.trajectory.tail.probs(sv.indel.nonsense.muts, aerobic.genes$Gene)
anaerobic.sv.indel.nonsense.pvals <- calculate.trajectory.tail.probs(sv.indel.nonsense.muts, anaerobic.genes$Gene)

######
## plot structural variation, indels, nonsense mutations.

## aerobic in black
c.aerobic.sv.indel.nonsense.muts <- sv.indel.nonsense.muts %>%
    filter(aerobic==TRUE) %>%
    calc.cumulative.muts()

## anaerobic in red
c.anaerobic.sv.indel.nonsense.muts <- sv.indel.nonsense.muts %>%
    filter(anaerobic==TRUE) %>%
    calc.cumulative.muts()

## make and save plot.
aerobic.anaerobic.purifying.plot <- plot.base.layer(sv.indel.nonsense.muts) %>%
    add.cumulative.mut.layer(c.aerobic.sv.indel.nonsense.muts, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.sv.indel.nonsense.muts, my.color="red")
ggsave(filename="../results/figures/SVIndelNonsense.pdf", plot=aerobic.anaerobic.purifying.plot)

######
## examine dN.
gene.dN.mutation.data <- gene.mutation.data %>%
filter(Annotation=='missense')
## calculate p-values for dN.
aerobic.dN.pvals <- calculate.trajectory.tail.probs(gene.dN.mutation.data, aerobic.genes$Gene)
anaerobic.dN.pvals <- calculate.trajectory.tail.probs(gene.dN.mutation.data, anaerobic.genes$Gene)

## plot just dN.
c.aerobic.dN.mutations <- gene.dN.mutation.data %>%
    filter(aerobic==TRUE) %>%
    calc.cumulative.muts()
c.anaerobic.dN.mutations <- gene.dN.mutation.data %>%
    filter(anaerobic==TRUE) %>%
    calc.cumulative.muts()

aerobic.anaerobic.dN.plot <- plot.base.layer(gene.dN.mutation.data) %>%
        add.cumulative.mut.layer(c.aerobic.dN.mutations, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.dN.mutations, my.color="red")    
ggsave(filename="../results/figures/dN.pdf", plot=aerobic.anaerobic.dN.plot)


######
## examine dS over the genome.
gene.dS.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='synonymous')
## calculate p-values for dS.
aerobic.dS.pvals <- calculate.trajectory.tail.probs(gene.dS.mutation.data, aerobic.genes$Gene)
anaerobic.dS.pvals <- calculate.trajectory.tail.probs(gene.dS.mutation.data, anaerobic.genes$Gene)

## plot dS.
c.aerobic.dS.mutations <- gene.dS.mutation.data %>%
    filter(aerobic==TRUE) %>%
    calc.cumulative.muts()

c.anaerobic.dS.mutations <- gene.dS.mutation.data %>%
    filter(anaerobic==TRUE) %>%
    calc.cumulative.muts()

aerobic.anaerobic.dS.plot <- plot.base.layer(gene.dS.mutation.data) %>%
    add.cumulative.mut.layer(c.aerobic.dS.mutations, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.dS.mutations, my.color="red") 
ggsave(filename="../results/figures/dS.pdf", plot=aerobic.anaerobic.dS.plot)

######
## nonsense mutations.
gene.nonsense.mutation.data <- gene.mutation.data %>%
    filter(Annotation=='nonsense')
## calculate p-values for nonsense.
aerobic.nonsense.pvals <- calculate.trajectory.tail.probs(gene.nonsense.mutation.data, aerobic.genes$Gene)
anaerobic.nonsense.pvals <- calculate.trajectory.tail.probs(gene.nonsense.mutation.data, anaerobic.genes$Gene)

## plot just nonsense mutations.
c.aerobic.nonsense.mutations <- gene.nonsense.mutation.data %>%
    filter(aerobic==TRUE) %>%
    calc.cumulative.muts()
c.anaerobic.nonsense.mutations <- gene.nonsense.mutation.data %>%
    filter(anaerobic==TRUE) %>%
    calc.cumulative.muts()

aerobic.anaerobic.nonsense.plot <- plot.base.layer(gene.nonsense.mutation.data) %>%
    add.cumulative.mut.layer(c.aerobic.nonsense.mutations, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.nonsense.mutations, my.color="red") 
ggsave(filename="../results/figures/nonsense.pdf", plot=aerobic.anaerobic.nonsense.plot)

######
## summary of p-value results.

p1 <- aerobic.allmut.pvals %>% mutate(mut.class='all') %>% mutate(gene.class='aerobic')
p2 <- anaerobic.allmut.pvals %>% mutate(mut.class='all') %>% mutate(gene.class='anaerobic')

p3 <- aerobic.sv.indel.nonsense.pvals %>% mutate(mut.class='sv.indel.nonsense') %>% mutate(gene.class='aerobic')
p4 <- anaerobic.sv.indel.nonsense.pvals %>% mutate(mut.class='sv.indel.nonsense') %>% mutate(gene.class='anaerobic')

p5 <- aerobic.dN.pvals %>% mutate(mut.class='dN') %>% mutate(gene.class='aerobic')
p6 <- anaerobic.dN.pvals %>% mutate(mut.class='dN') %>% mutate(gene.class='anaerobic')

p7 <- aerobic.dS.pvals %>% mutate(mut.class='dS') %>% mutate(gene.class='aerobic')
p8 <- anaerobic.dS.pvals %>% mutate(mut.class='dS') %>% mutate(gene.class='anaerobic')

p9 <- aerobic.nonsense.pvals %>% mutate(mut.class='nonsense') %>% mutate(gene.class='aerobic')
p10 <- anaerobic.nonsense.pvals %>% mutate(mut.class='nonsense') %>% mutate(gene.class='anaerobic')

pval.summary <- rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10) %>%
    mutate(hypermutator=ifelse((Population %in% hypermutator.pops),TRUE,FALSE))
    
write.csv(pval.summary,"../results/trajectory-pvalues.csv")

######
## let's look at noncoding mutations.
## WARNING: we have to be careful about the appropriate normalization for noncoding mutations.
## I am assuming that these are within the ORF-- but I haven't explicitly checked.
## Hold off on report this result for now (maybe cut entirely.)

## divide by gene length since we've filtered out intergenic mutations
## from gene.mutation.data.
## Still not sure if this is the right normalization.

gene.noncoding.mutation.data <- filter(gene.mutation.data,Annotation=='noncoding')

c.aerobic.noncoding.mutations <- gene.noncoding.mutation.data %>%
    filter(aerobic==TRUE) %>%
    calc.cumulative.muts(normalization.constant=NA)

c.anaerobic.noncoding.mutations <- gene.noncoding.mutation.data %>%
    filter(anaerobic==TRUE) %>%
    calc.cumulative.muts(normalization.constant=NA)

aerobic.anaerobic.noncoding.plot <- plot.base.layer(gene.noncoding.mutation.data,
                                                    normalization.constant=NA) %>%
    add.cumulative.mut.layer(c.aerobic.noncoding.mutations, my.color="black") %>%
    add.cumulative.mut.layer(c.anaerobic.noncoding.mutations, my.color="red")
aerobic.anaerobic.noncoding.plot
ggsave(filename="../results/figures/noncoding.pdf", plot=aerobic.anaerobic.noncoding.plot)

########################################################################################
## LOOK AT PARALLEL EVOLUTION IN NON-CODING MUTATIONS!
## LOOK AT mntH and hupA
aerobic.noncoding.summary <- gene.noncoding.mutation.data %>%
    filter(aerobic==TRUE) %>% group_by(Gene,Position) %>%
    summarize(count=n()) %>% arrange(desc(count,Position))

aerobic.noncoding.summary2 <- gene.noncoding.mutation.data %>%
    filter(aerobic==TRUE) %>% group_by(Gene) %>%
    summarize(count=n()) %>% filter(count>1) %>% arrange(desc(count,Position))

aerobic.noncoding.summary3 <- filter(aerobic.noncoding.summary, Gene %in% aerobic.noncoding.summary2$Gene)

anaerobic.noncoding.summary <- gene.noncoding.mutation.data %>%
    filter(anaerobic==TRUE) %>% group_by(Gene,Position) %>% summarize(count=n()) %>% arrange(desc(count,Position))

anaerobic.noncoding.summary2 <- gene.noncoding.mutation.data %>%
    filter(anaerobic==TRUE) %>% group_by(Gene) %>%
    summarize(count=n()) %>% filter(count > 1) %>% arrange(desc(count,Position))

anaerobic.noncoding.summary3 <- filter(anaerobic.noncoding.summary, Gene %in% anaerobic.noncoding.summary2$Gene)

## let's look at all noncoding mutations-- not just those associated with genes.
noncoding.mutation.data <- mutation.data %>% filter(Annotation=='noncoding')

noncoding.summary1 <- noncoding.mutation.data %>% group_by(Gene,Position) %>%
    summarize(count=n()) %>% arrange(desc(count,Position))

noncoding.summary2 <- noncoding.mutation.data %>% group_by(Gene) %>%
    summarize(count=n()) %>% filter(count > 3) %>% arrange(desc(count,Position))

noncoding.summary3 <- filter(noncoding.summary1, Gene %in% noncoding.summary2$Gene) 

noncoding.summary4 <- filter(noncoding.summary3,Gene != "intergenic")
noncoding.summary4

mntH.muts <- filter(mutation.data,Gene == 'mntH')
