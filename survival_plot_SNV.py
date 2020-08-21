import matplotlib.pyplot as plt
import numpy as np

def parse_data(filename):
    file = open(filename,"r")
    file.readline() # header
    #Cohort, Sample_t0, Sample_t1, Species, NumTestedSNPs, SNPs, NumTestedGenes
    snp_data_hmp = []
    snp_data_twins = []
    snp_data_youngtwins = []
    
    for line in file:
        items = line.split(",")
        
        cohort = items[0].strip();
        num_snp_changes = int(items[5].strip());

        
        if cohort == "hmp":
            snp_data_hmp.append(num_snp_changes)
        if cohort == "twins":
            snp_data_twins.append(num_snp_changes)
        if cohort == "young_twins":
            snp_data_youngtwins.append(num_snp_changes)    
        
        
    return snp_data_hmp, snp_data_twins, snp_data_youngtwins

def create_figure(snp_data_hmp, snp_data_twins, snp_data_youngtwins):   
    modification_max = 20;
    replacement_min = 500
    
    fig, ax = plt.subplots(figsize=(7, 4))
      
    n_bins = 500000 
    
    #shaded regions
    ax.axvspan(0, 1, alpha=0.5, color='#818281')
    ax.axvspan(1, modification_max, alpha=0.5, color='#9fc9fc')
    ax.axvspan(replacement_min, 100000, alpha=0.5, color='#ffbe9c')
    
    #label shaded regions
    ax.text(
        10000, 1.3, 'putative\nreplacement', fontsize=13, ha='center',
        color='#ffbe9c', fontstyle='italic',
    )
    ax.text(
        5, 1.3, 'putative\nmodification', fontsize=13, ha='center',
        color='#9fc9fc', fontstyle='italic',
    )
        
        
    
    #plot data
    ax.hist(
        snp_data_hmp, n_bins, density=True, histtype='step',cumulative=-1, 
        label='within-host\n(hmp)', color = '#0052a3'
    )
    ax.hist(
        snp_data_twins, n_bins, density=True, histtype='step',cumulative=-1, 
        label='between-host\n(adult twins)', color = 'purple'
    )
    ax.hist(
        snp_data_youngtwins, n_bins, density=True, histtype='step',cumulative=-1, 
        label='between-host\n(young twins)', color = 'green'
    )    
    
    ax.legend(bbox_to_anchor=(1,1)) #change to place legend outside of plot
    
    plt.subplots_adjust(right=0.7)
    
    #axis labels, range
    ax.set_xlabel("# SNV changes")
    ax.set_ylabel("Fraction comparisons $≥ n$")
    ax.set_yscale('log')    
    ax.set_xscale('log') 
    ax.set_xlim(.6, 105000)
    ax.set_ylim(1/800, 1.1) 
    
    #ticks
    ax.tick_params(which = "both", direction = "in")
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=3)
    
    
    #plot outlines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)    
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    
    plt.show() 
    plt.savefig("figure1")
    
    
if __name__=='__main__': 
    snp_data_hmp, snp_data_twins, snp_data_youngtwins = parse_data("within_host_changes.csv")
    create_figure(snp_data_hmp, snp_data_twins, snp_data_youngtwins)
    
    