import numpy as np
import matplotlib.pyplot as plt
import random
import csv
from matplotlib.ticker import MaxNLocator



def parse_data(filename):
    
    file = open(filename,"r")
    file.readline() # header
    
    
    data = [] #Cohort, Sample_t0, Species,SNPs
  
    for line in file:
        items = line.split(",")
        
        cohort = items[0].strip();
        sample_t0 = items[1].strip();
        sample_t1 = items[2].strip();
        
        species = items[3].strip();        
        num_snp_changes = int(items[5].strip());
        
        data.append((cohort, species, (sample_t0,sample_t1), num_snp_changes))
        
    return data
        

def sort_data(data):
    THRESHOLD = 500
    
    hmp_dict = {}
    twins_dict = {}
    y_twins_dict = {}   
    #for each dictionary: the key is Sample_t0 and the value is a list containing information on the sample
    #the list (by host): num species, num species with rep OR mod, num species with mod, num species with rep
    
    
    for species in data:
        cohort = species[0]
        sample_t0 = species[2]
        num_snp_changes = species[3]
        
        if cohort == "hmp":
            d = hmp_dict
        elif cohort == "twins":
            d = twins_dict
        elif cohort == "young_twins":
            d = y_twins_dict
        else:
            #THROW ERROR
            print(cohort)
            print("not in 3 data sets")
        
        
        if sample_t0 not in d:
            d[sample_t0] = [1,0,0,0]
        else:
            d[sample_t0][0]+=1
            
        if num_snp_changes!=0:
            d[sample_t0][1]+=1
            if num_snp_changes > THRESHOLD:
                d[sample_t0][3]+=1
            else:
                d[sample_t0][2]+=1        
                            
    return hmp_dict, twins_dict, y_twins_dict
    

def create_figure(hmp, twins, y_twins, labels):
    
    hmp_rep_list = []
    hmp_mod_list = []
    twins_rep_list = []
    twins_mod_list = []
    y_twins_rep_list = []
    y_twins_mod_list = []
    total_hmp_list = []
    total_twins_list = []
    total_y_twins_list = []
    hmp_or_list = []
    twins_or_list = []
    y_twins_or_list = []    
    
    hmp_ratio = []
    twins_ratio = []
    y_twins_ratio = []
    
    for id in hmp:
        total_hmp_list.append(hmp[id][0])
        hmp_or_list.append(hmp[id][1])
        hmp_mod_list.append(hmp[id][2])
        hmp_rep_list.append(hmp[id][3])
        hmp_ratio.append(hmp[id][1]/hmp[id][0])
        
    for id in twins:
        total_twins_list.append(twins[id][0])
        twins_or_list.append(twins[id][1])
        twins_mod_list.append(twins[id][2])
        twins_rep_list.append(twins[id][3]) 
        twins_ratio.append(twins[id][1]/twins[id][0])
        
    for id in y_twins:
        total_y_twins_list.append(y_twins[id][0])
        y_twins_or_list.append(y_twins[id][1])
        y_twins_mod_list.append(y_twins[id][2])
        y_twins_rep_list.append(y_twins[id][3])
        y_twins_ratio.append(y_twins[id][1]/y_twins[id][0])
    
    #count = 0
    #for id in hmp:
        #if (hmp[id][2] > 0) and (hmp[id][3] > 0):
            #count += 1
    #print("number of hosts with modification and replacement events : " + str(count))
    

    BINS = [0, 1, 2, 3, 4, 5, 6]
    
    fig, axs = plt.subplots(3, 2, figsize=(10, 7), sharey='row')
    fig.delaxes(axs[2,1])
    
    axs[0][0].set_yscale('log')
    axs[1][0].set_yscale('log') #maybe
         
    
    colors = ['#5482ff', '#b175ff', '#79ba92']
    
    #modification events
    axs[0][0].hist(
        [hmp_mod_list, twins_mod_list, y_twins_mod_list], bins=BINS, histtype='bar', 
        color=colors, label=labels, align='left'
    )  
    axs[0][0].set_xlabel("# of species with modification events within host")
    axs[0][0].set_ylabel("# hosts") 
    
    
    #replacement events
    axs[0][1].hist(
        [hmp_rep_list, twins_rep_list, y_twins_rep_list], bins=BINS, histtype='bar', 
        color=colors, label=labels, align='left'
    )  
    axs[0][1].set_xlabel("   # of species with replacement events within host")
    axs[0][0].set_ylim(.9,200)
    
    axs[1][0].hist(
        [hmp_or_list, twins_or_list, y_twins_or_list], bins=BINS, histtype='bar', 
        color=colors, label=labels, align='left'
    )
    axs[1][0].set_xlabel("# of species with replacement or modifications within host")
    axs[1][0].set_ylabel("# hosts")
    
    axs[1][0].set_ylim(.9,200)
    
    
    
    axs[1][1].hist(
        [total_hmp_list, total_twins_list, total_y_twins_list], bins = 12, histtype='bar', 
        color=colors, label=labels, align='left'
    )
    axs[1][1].set_xlabel("# of species in host")
    axs[1][1].xaxis.set_ticks([0,2,4,6,8,10,12])
    
    axs[2][0].set_yscale('log')
    axs[2][0].set_ylim(1/200, 1.5)
    axs[2][0].set_xlabel("fraction of species in host with replacement or modification events")
    axs[2][0].set_ylabel("Fraction comparisons $≥ n$")

    
    axs[2][0].hist(
        [hmp_ratio,twins_ratio,y_twins_ratio], 500000, density=True, histtype='step',cumulative=-1, 
        label=labels, color = colors
    )    

    axs[0][0].text(-0.02, 1.1, 'A', transform=axs[0][0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    axs[0][1].text(-0.02, 1.1, 'B', transform=axs[0][1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    axs[1][0].text(-0.02, 1.1, 'C', transform=axs[1][0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    axs[1][1].text(-0.02, 1.1, 'D', transform=axs[1][1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    axs[2][0].text(-0.02, 1.1, 'E', transform=axs[2][0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
        
    
    
    axs[0][1].legend(bbox_to_anchor=(1.5,1)) #change to place legend outside of plot
    plt.subplots_adjust(right=0.8, top = .95, left = 0.1, bottom = .1, wspace = 0.15, hspace = .4 )
    
    
    

    plt.show() 
    plt.savefig("figures/observed_data")    

def generate_fake_data_set(filename, species_permutations):
    
    
    hmp_data = []
    snv_list = []
    
    file = open(filename,"r")
    file.readline()     
    
    if species_permutations:
        species_dict = {}
    
    for line in file:
        items = line.split(",")
        
        cohort = items[0].strip();
        sample_t0 = items[1].strip();
        species = items[3].strip();
        num_snp_changes = int(items[5].strip());
        
        if cohort == "hmp":
            l = [cohort, species, sample_t0]
            hmp_data.append(l)
            
            if species_permutations:
                if species not in species_dict:
                    species_dict[species] = []
                species_dict[species].append(num_snp_changes)
                    
            else:
                snv_list.append(num_snp_changes)
                
    
    if species_permutations:
        for s in species_dict:
            random.shuffle(species_dict[s])        

        #for each sample in the hmp data
        for sample in hmp_data:
            #append to the sample a SNP value that corresponds to the 
            species_name = sample[1]
            sample.append(species_dict[species_name][0])
            species_dict[species_name].pop(0)
                    
        
    else:
        random.shuffle(snv_list)
        
        i = 0
        for species in hmp_data:
            species.append(snv_list[i])
            i += 1
            
    return hmp_data

def calculate_stats(null_data):

    rep_list = []
    mod_list = []
    total_list = []
    or_list = []
    
    for id in null_data:
        total_list.append(null_data[id][0])
        or_list.append(null_data[id][1])
        mod_list.append(null_data[id][2])
        rep_list.append(null_data[id][3])
            
    events_2 = 0
    events_3 = 0
    mod_2 = 0
    mod_3 = 0
    rep_2 = 0
    rep_3 = 0
    
    for num in or_list:
        if num >= 2:
            events_2 += 1 
            if num >= 3:
                events_3 += 1
    
    for num in mod_list:
        if num >= 2:
            mod_2 += 1
            if num>=3:
                mod_3 += 1
                
    for num in rep_list:
        if num >= 2:
            rep_2 += 1
            if num>= 3:
                rep_3 += 1    

    stats = [events_2, events_3, mod_2, mod_3, rep_2, rep_3]
    
    for i in range(0,len(stats)):
        stats[i] = float(stats[i])
                    
    return stats
  
    
def save_data(null_stats, null_filename):
    with open(null_filename, mode='a') as file:
        file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        file_writer.writerow(null_stats)


def run_real(filename):
    data_real = parse_data(filename)
    hmp, twins, y_twins = sort_data(data_real)
    create_figure(hmp, twins, y_twins, ["hmp", "twins", "young_twins"])        


def make_null_data(experimental_data_filename, null_data_filename, num_data_sets, species = True):
    
    random.seed(0)
    
    with open(null_data_filename, mode='w') as file:
        file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)        
        file_writer.writerow(["events>=2", "events>=3", "mod>=2", "mod>=3", "rep>=2", "rep>=3"])
        
    
    for i in range(0,num_data_sets):
        null_data = generate_fake_data_set(experimental_data_filename, species)
        sorted_null_data, a, b, = sort_data(null_data)
        null_stats = calculate_stats(sorted_null_data)
        save_data(null_stats, null_data_filename)


def display_null_distributions(experimental_data_filename, null_data_filename):
    data_real = parse_data(experimental_data_filename)
    hmp, twins, y_twins = sort_data(data_real)
    
    observed_vals = calculate_stats(hmp)
    
    file = open(null_data_filename,"r")
    file.readline() # header
    
        
    events_2 = []
    events_3 = []
    mod_2 = []
    mod_3 = []
    rep_2 = []
    rep_3 = []
    
    for line in file:
        items = line.split(",")
        events_2.append(float(items[0].strip()))
        events_3.append(float(items[1].strip()))
        mod_2.append(float(items[2].strip()))
        mod_3.append(float(items[3].strip()))
        rep_2.append(float(items[4].strip()))
        rep_3.append(float(items[5].strip()))
        
        
    #calculate p-values
    p_values = []
    stats_lists = [events_2, events_3, mod_2, mod_3, rep_2, rep_3]
    
    for i in range(0, len(stats_lists)):
        count = 0        
        for data_set in stats_lists[i]:
            if observed_vals[i] <= float(data_set):
                count += 1
        
        p_values.append(((count+1)/(len(stats_lists[i])+1)))
   
    fig, axs = plt.subplots(3, 2, figsize=(10, 7))
    
    BINS = 100
    titles = ["# hosts with 2 or more modification or replacement events", 
              "# of hosts with 3 or more modification or replacement events",
              "# of hosts with 2 or more modification events", 
              "# of hosts with 3 or more modification events",
              "# of hosts with 2 or more replacement events",
              "# of hosts with 3 or more replacement events"
              ]
    
    data_list = [events_2, events_3, mod_2, mod_3, rep_2, rep_3]
    stats_lists = [events_2, events_3, mod_2, mod_3, rep_2, rep_3]

     
    p_line_start_x = [10, 0, 5,0,0,0]
    lab = ['A', 'B', 'C', 'D', 'E', 'F']
    i = 0
    for r in range(0,3):
        for c in range(0,2):
            axs[r][c].hist(data_list[i],density=True, histtype='step',cumulative=-1, bins=BINS,color='grey', label = "simulated hmp data", linewidth=2)
            
            axs[r][c].set_yscale('log')
            axs[r][c].set_ylim(1/6000,1.5)
            axs[r][c].set_ylabel("Fraction comparisons $≥ n$")
            axs[r][c].set_xlabel(titles[i])
            
            axs[r][c].xaxis.set_major_locator(MaxNLocator(integer=True))
           
            axs[r][c].plot([p_line_start_x[i], observed_vals[i]], [p_values[i], p_values[i]], color = 'red', ls = '--', linewidth=1, label = "observed value, p = " + str(p_values[i]))
            axs[r][c].plot([observed_vals[i], observed_vals[i]], [1/6000, p_values[i]], color = 'red', ls = '--', linewidth=1)    
            axs[r][c].scatter(observed_vals[i], p_values[i], s = 30, color = 'red')
            
            axs[r][c].text(0, 1.15, lab[i], transform=axs[r][c].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
                        
            
            axs[r][c].legend() 
            
            
            i += 1
        
    plt.subplots_adjust(right=0.95, top = .95, left = 0.07, bottom = .08, wspace = 0.32, hspace = .53 )
    
    plt.show()
    plt.savefig("figures/null_distributions_species") 
    

if __name__=='__main__': 
    
    run_real("within_host_changes.csv")
    make_null_data("within_host_changes.csv","null_data_species.csv", 9999, False)
    display_null_distributions("within_host_changes.csv","null_data_species.csv")
    
   