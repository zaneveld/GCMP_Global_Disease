import argparse
import pathlib
from os.path import join,exists
from os import listdir
import subprocess
from pandas import DataFrame
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
from numpy import array
import numpy as np
import matplotlib.pyplot as plt
from  numpy import array


def make_argument_parser():
    """Return an argument parser for command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--trait_table_fp", type=pathlib.Path,
                    required = True,
                    help="input trait table filepath")
    parser.add_argument("-t","--tree_fp",type=pathlib.Path, \
                    required = True,\
                    help="phylogenetic tree. Tips must match taxa in trait table.")
    parser.add_argument("-o","--output_results_dir",required=True,\
                    default=None,type=pathlib.Path,\
                    help="top level results directory for output. \
                    Subdirectory will be created within based on the\
                    analysis name. Default: ../output")
    parser.add_argument('-c','--host_trait_column',type=str,required=True,\
                    help = "Name of column that holds the host trait \
                    that you want to regress against microbial relative abundance")
    parser.add_argument("--verbose", help="print verbose output",
                    action="store_true",required=False)
    parser.add_argument("-l","--analysis_label",default="taxon_bubble_plot",
                    required=False)
    parser.add_argument("--test_run",default=False,\
                    action="store_true",required=False)
    return parser


def main():
    parser = make_argument_parser()
    args = parser.parse_args()
    print(args)
    trait_table_fp = args.trait_table_fp
    results_dir = args.output_results_dir or join("..","output")
    tree_fp = args.tree_fp
    analysis_label = args.analysis_label
    host_trait = args.host_trait_column
    analysis_output_dir = join(results_dir, analysis_label)

    if args.verbose:
        print("Trait table:",trait_table_fp)
        print("Tree:",tree_fp)
        print("Results dir:",results_dir)
        print("Host trait:",host_trait)
        print("Analysis output dir:", analysis_output_dir)


    dominant_genera,prevalence = retrieve_dominant_genera(trait_table_fp,
      results_dir,analysis_label="growth_vs_microbe_genera")

    if args.test_run:
        print("Truncating dominant genera to 5 for test run!")
        dominant_genera = dominant_genera[0:5]

    short_names_for_dominant_genera = [d.split("__")[3]+"__"+d.split("__")[-1] for d in dominant_genera]
    print(f"Dominant genera: {short_names_for_dominant_genera}")

    #Run PIC and PGLS
    pic_results_df,pgls_results_df,bad_genera = run_trait_vs_taxonomy_PGLS(trait_table_fp,tree_fp,\
      target_taxa=dominant_genera,analysis_label = analysis_label,\
       trait = host_trait,results_dir = results_dir,prevalence = prevalence)

    #Save results dataframe
    pic_results_df.to_csv(join(analysis_output_dir,\
      f"PIC_results_summary_{analysis_label}.tsv"),sep="\t")

    pgls_results_df.to_csv(join(analysis_output_dir,\
      f"PGLS_results_summary_{analysis_label}.tsv"),sep="\t")

    #Parse PGLS results to get graph data
    x_values,y_values,colors,sizes,markers,labels,compartments,\
    x_position_mapping,y_position_mapping,color_mapping =\
      get_taxon_bubble_plot_data(analysis_output_dir,analysis_label)

    #Generate and save graph
    make_taxonomy_bubble_plot(analysis_output_dir,analysis_label,\
      x_values,y_values,colors,sizes,markers,labels,compartments,\
      x_position_mapping,y_position_mapping,color_mapping)





def run_trait_vs_taxonomy_PGLS(trait_table_fp,tree_fp,target_taxa,\
  analysis_label,trait,results_dir,prevalence = None):
    bad_genera = []
    pic_results_df = DataFrame({},columns = ["analysis_label","pic_x_trait","pic_y_trait","R2","p","sig_marker","FDR_q","slope","pic_filter_column","pic_filter_value","results_dir","slope_std_error","T_stat"])
    pgls_results_df = DataFrame({},columns = ["analysis_label","x_trait","y_trait","R2","p","FDR_q","N_unique_samples","slope","model_name","best_model","AIC","AICc","delta_AICc","filter_column","filter_value","results_dir","x_trait_slope_95CI"])
    if not prevalence:
        prevalence = ['Not calculated']*len(target_taxa)
    dominant_genera = target_taxa
    pic_y_trait = trait
    analysis_output_dir = join(results_dir, analysis_label)
    for i,dominant_genus in enumerate(dominant_genera):
        pic_tree = tree_fp
        pic_trait_table = trait_table_fp
        pic_x_trait = dominant_genus
        #pic_x_trait = pic_x_trait.replace(";","__").replace(":","_")


        pic_filter_column = 'None'
        pic_filter_value = 'None'

        pic_suffix = ''
        try:
            pic_result = phylogenetic_independent_contrasts(pic_trait_table,\
             pic_tree, pic_x_trait,pic_y_trait,pic_filter_column,\
            pic_filter_value,analysis_output_dir,verbose = False)
            pic_result["analysis_label"] = analysis_label
            pic_results_df = pic_results_df.append(pic_result,ignore_index=True)
            pic_results_df = add_FDR(pic_results_df)
        except TypeError:
            print("Skipping genus:",dominant_genus, "due to PIC error (typically no variance / absent in corals that map to the tree)")
            bad_genera.append(dominant_genus)
            print("Skipped genera:",bad_genera)
            print(f"Done! Skippped {len(bad_genera)} for which PIC or PGLS could not be calculated")
            continue

        try:
            pgls_result = pgls(pic_trait_table,pic_tree,pic_x_trait,pic_y_trait,pic_filter_column,\
              pic_filter_value,analysis_output_dir,verbose = False)
            pgls_result["analysis_label"] = analysis_label
            pgls_result["N_unique_samples"] = prevalence[i] #<- calculated in previous cell
            pgls_results_df = pgls_results_df.append(pgls_result,ignore_index=True)
            pgls_results_df = add_FDR(pgls_results_df,best_model_only = True)
        except TypeError:
            print("Skipping genus:",dominant_genus, "due to PGLS error (typically no variance / absent in corals that map to the tree)")
            bad_genera.append(dominant_genus)
            print("Skipped genera:",bad_genera)
            print(f"Done! Skippped {len(bad_genera)} for which PIC or PGLS could not be calculated")
            continue




    return pic_results_df,pgls_results_df,bad_genera

def phylogenetic_independent_contrasts(pic_trait_table,pic_tree,pic_x_trait,\
    pic_y_trait,pic_filter_column,pic_filter_value,\
    output_dir="../PIC_results/",pic_suffix='',verbose=True):
    """Run the phylomorphospace script"""
    #Build up the command we want to run
    pic_cmd = f"Rscript phylomorphospace_r14.r {pic_trait_table} {pic_tree} {pic_x_trait} {pic_y_trait} {pic_filter_column} {pic_filter_value} {output_dir} {pic_suffix}"
    print("Running PIC command:",pic_cmd)

    try:
        pic_output = subprocess.check_output(pic_cmd.split(),stderr=subprocess.STDOUT)
        pic_output = str(pic_output)
    except subprocess.CalledProcessError as exc:
        print(exc.output)
        return exc.output

    result_lines = pic_output.split("\\n")
    if verbose:
        for line in result_lines:
            print(line)

    results = parse_pic_result_lines(result_lines)
    results['trait_table'] = pic_trait_table
    results['tree'] = pic_tree
    results['pic_x_trait'] = pic_x_trait
    results['pic_y_trait'] = pic_y_trait
    results['pic_filter_column'] = pic_filter_column
    results['pic_filter_value'] = pic_filter_value
    results['pic_suffix'] = pic_suffix

    return results

def pgls(trait_table,tree,x_trait,y_trait,filter_column,filter_value,output_dir="../PIC_results/",verbose=True):
    """Run the PGLS script"""

    #Build up the command we want to run
    pgls_cmd = f"Rscript pgls.r {trait_table} {tree} {x_trait} {y_trait} {filter_column} {filter_value} {output_dir}"
    print("Running PGLS script as follows:",pgls_cmd)

    try:
        pgls_output = subprocess.check_output(pgls_cmd.split(),stderr=subprocess.STDOUT)
        pgls_output = str(pgls_output)
    except subprocess.CalledProcessError as exc:
        print(exc.output)
        return exc.output

    result_lines = pgls_output.split("\\n")
    if verbose:
        for line in result_lines:
            print(line)

    output_filepath = None
    for line in result_lines:
        if "Writing PGLS results file" in line:
            output_filepath = line.split(":")[1].rstrip('"')
            print("Output filepath:",output_filepath)
    if not output_filepath:
        raise ValueError("Can't find output filepath in results")

    results = parse_pgls_results(output_filepath)

    #Add compartment column
    #NOTE: the PGLS script is run once per compartment,
    #so at this stage we can safely assume that there
    #is only 1 compartment represented in the results
    compartments = ['all','mucus','tissue','skeleton']
    results['compartment'] = 'None'
    for possible_compartment in compartments:
        x_trait_values = results['x_trait']
        for i,x_trait_value in enumerate(x_trait_values):
            if x_trait_value.endswith(possible_compartment) or\
              x_trait_value.startswith(possible_compartment):
                results.loc[i,"compartment"] = possible_compartment

    #results = parse_pic_result_lines(result_lines)
    results['trait_table'] = trait_table
    results['tree'] = tree
    results['x_trait'] = x_trait
    results['y_trait'] = y_trait
    results['filter_column'] = filter_column
    results['filter_value'] = filter_value
    #results['pic_suffix'] = pic_suffix


    return results

def parse_pgls_results(pgls_results_filepath):
    """Parse PGLS results and return the result
    """
    df = pd.read_csv(pgls_results_filepath,sep="\t")
    df = add_delta_AICc(df)
    df['results_dir'] = pgls_results_filepath
    return df

def add_delta_AICc(df,aicc_column_name='AICc'):
    "Add a DataFrame column with delta AICc values"
    best_aic = min(df[aicc_column_name])
    df[f'delta_{aicc_column_name}'] = df[aicc_column_name] - best_aic
    return df


def parse_pic_result_lines(lines):
    results = {}
    for line in lines:
        if line.startswith("pic.X"):
            fields = line.split()[1:]
            if len(fields)==4:
                slope,std_error,T,p =  fields
                sig_marker = 'n.s'
            else:
                slope,std_error,T,p,sig_marker = fields


            results['slope'] = slope
            results['slope_std_error'] = std_error
            results['T_stat'] = T
            results['p'] = p
            results['sig_marker'] = sig_marker

        if line.startswith('[1] \"Outputting results to: '):
            results['results_dir'] = line.split(":")[1].rstrip(",")
        if line.startswith("Multiple R-squared:"):
            R2 = float(line.split(":")[1].split(",")[0])
            results['R2'] = R2

    print("R2:",results['R2'])
    print("p:",p)
    return results



def add_FDR(df,p_value_column_name = "p",best_model_only = False):
    """Add FDR q values to a DataFrame in column FDR_q

    df -- pandas DataFrame
    p_value_column_name -- the column of the DataFrame holding p p_values
    best_model_only -- if True, use only rows where a 'best_model' columns
      is set to True.
    """
    if best_model_only:
        #ignore high (worse) AICc models
        df = df.sort_values('best_model',ascending=False)
        df = df.reset_index(drop=True)
        best_models = df[df['best_model'] == True]
        p_values = list(best_models[p_value_column_name])
    else:
        p_values = list(df[p_value_column_name])

    p_values = array(list(map(float,[p.strip("<") if type(p) is str else p for p in p_values])))
    rejected,fdr_values = fdrcorrection(p_values,alpha=0.05,method='indep',is_sorted=False)


    if best_model_only:
        df['FDR_q'] = 'NA (q values only calculated for best models by AIC)'

        for i,q in enumerate(fdr_values):
            df.loc[i,"FDR_q"] = q

    else:
        df['FDR_q'] = fdr_values
    return df

def retrieve_dominant_genera(trait_table,results_dir,analysis_label="growth_vs_microbe_genera"):

    pic_results_df = DataFrame({},columns = ["analysis_label","pic_x_trait","pic_y_trait","R2","p","sig_marker","FDR_q","slope","pic_filter_column","pic_filter_value","results_dir","slope_std_error","T_stat"])
    pgls_results_df = DataFrame({},columns = ["analysis_label","x_trait","y_trait","R2","p","FDR_q","N_unique_samples","slope","model_name","best_model","AIC","AICc","delta_AICc","filter_column","filter_value","results_dir","x_trait_slope_95CI"])

    analysis_output_dir = join(results_dir,"PIC_results",f"{analysis_label}")

    compartments = ["all","mucus","tissue","skeleton"]

    #Find taxon names
    pic_trait_table = pd.read_csv(trait_table,sep="\t")
    dominant_genera = []

    prevalence = []

    for compartment in compartments:
        for column in list(pic_trait_table.columns):
            if  ("D_5" in column) & (compartment in column) & ("D_6" not in column):
                #Filter to only microbes present in at least 3 genera
                data = list(pic_trait_table[column])
                n_above_zero = 0
                for d in data:
                    if float(d) > 0.0:
                        n_above_zero += 1
                #print(f"{column}: present in {n_above_zero} genera")

                if n_above_zero <= 3:
                    #print("Skipping ... too few observations")
                    continue

                dominant_genera.append(str(column))
                prevalence.append(str(n_above_zero))

    #Save only unique genera (no duplicates)
    dominant_genera = [g for g in dominant_genera if str(g) != 'nan']
    dominant_genera = list(set(dominant_genera))
    return dominant_genera,prevalence

def make_taxonomy_bubble_plot(analysis_output_dir,analysis_label,\
  x_values,y_values,colors,sizes,markers,labels,compartments,\
  x_position_mapping,y_position_mapping,color_mapping,ytick_scalar = 0.01):


    fig,ax = plt.subplots(figsize=(13,4))

    print(markers)
    print(x_values)
    print(y_values)
    print(colors)
    print(sizes)

    #Plot positive values
    for marker in ["^","v"]:
        print(markers)
        print(len(labels))
        print(labels)
        positive_x_values = [x for i,x in enumerate(x_values) if markers[i] == marker]
        positive_y_values = [y for i,y in enumerate(y_values) if markers[i] == marker]
        positive_sizes = [s for i,s in enumerate(sizes) if markers[i] == marker]
        positive_colors = [c for i,c in enumerate(colors) if markers[i] == marker]
        #positive_labels = [l for i,l in enumerate(labels) if markers[i] == marker]

        #ax.scatter(positive_x_values,positive_y_values,s=1,facecolor='k',edgecolor='k',alpha=1.0)

        ax.scatter(positive_x_values,positive_y_values,s=array(positive_sizes)*1.5,\
                   c=positive_colors,edgecolors='k',alpha=1.0,marker=marker)

        #print(marker, positive_labels)

    for compartment in compartments:
        ax.axhline(y_position_mapping[compartment]*ytick_scalar,\
          color=color_mapping[compartment],alpha=0.10,linewidth=45,zorder=0)

    #Set the custom x and y-axis labels to make this a bubble plot
    ax.set_yticks([i*ytick_scalar*0.5 for i in range(len(compartments))])
    ax.set_yticklabels([c.capitalize() for c in reversed(compartments)],size='large')

    ax.set_xticks([i for i in range(len(labels))])

    def format_xtick_labels(full_taxon_strings):
        result = []
        for t in full_taxon_strings:
            fields = t.split("__")
            name = fields[-7]+"; "+fields[-1]
            name = name.replace("_"," ")
            result.append(name)
            #result.append(fields[-1])
        return result
    #here
    ax.set_xticklabels([l for l in format_xtick_labels(labels)],size='medium',\
                         rotation = 45, ha = "right",rotation_mode="anchor")


    xtick_labels = format_xtick_labels(labels)
    ax.set_axisbelow(True)
    unique_x_values = list(range(len(labels)))

    margin = 0.005
    ax.set_ylim(min(y_values)-margin,max(y_values)+margin)
    ax.set_ylabel('Compartment',size = 'xx-large')
    ax.set_xlabel('Taxon',size = 'xx-large')
    ax.grid(color='lightgrey',alpha=0.25)

    #Special formatting for some labels
    for label in (ax.get_xticklabels()):
        if "Endozoicomonas" in str(label):
            label.set_weight("heavy")

    for label in (ax.get_yticklabels()):
        if "All" in str(label):
            label.set_weight("heavy")
            label.set_color("k")

    ## Remove axis lines.
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    fig.savefig(join(analysis_output_dir,f"taxon_bubble_plot_{analysis_label}.pdf"),bbox_inches="tight")
    fig.savefig(join(analysis_output_dir,f"taxon_bubble_plot_{analysis_label}.png"),bbox_inches="tight")


def get_taxon_bubble_plot_data(analysis_output_dir,analysis_label):
    pgls_result_file = join(analysis_output_dir,\
      f"PGLS_results_summary_{analysis_label}.tsv")
    pgls_df = pd.read_csv(pgls_result_file,sep="\t")

    print(list(pgls_df['best_model']))
    pgls_df = pgls_df[pgls_df['best_model'] == True]


    #Assume that taxonomy was used as the x-axis trait
    unique_taxa = list(pgls_df['x_trait'].unique())

    compartments = ['all','mucus','tissue','skeleton']
    compartment_results_dict = {}
    for compartment in compartments:

        #Set up an inner dict to hold per-organism results
        compartment_results_dict[compartment] = {}

        #Filter df to just the current compartment
        curr_pgls_df = pgls_df[pgls_df['compartment'] == compartment]

        for unique_taxon in unique_taxa:
            #print(unique_taxon)
            taxon_pgls_df = curr_pgls_df[curr_pgls_df['x_trait'] == unique_taxon]

            if taxon_pgls_df.empty:
                continue

            taxon_pgls_df = taxon_pgls_df.sort_values('AICc')
            best_model = taxon_pgls_df.iloc[0,:]
            print(dict(best_model))
            Rsquared = best_model['R2']
            compartment_results_dict[compartment][unique_taxon.split("_",1)[1]] = dict(best_model)

    results = pd.DataFrame(compartment_results_dict)

    unique_taxa_names = list(set([t.split("_",1)[1] for t in unique_taxa]))
    print(unique_taxa_names)

    x_position_mapping = {t:i for i,t in enumerate(sorted(unique_taxa_names))}
    y_position_mapping = {'all':1.5,'mucus':1.0,'tissue':0.5,'skeleton':0}

    color_mapping = {'all':'black','mucus':'cyan','tissue':'orange','skeleton':'purple'}
    ytick_scalar = 0.01
    x_values = []
    y_values = []
    sizes    = []
    labels   = []
    colors   = []
    shapes   = []
    sig = []
    markers = []

    #Exclude taxa not found in at least min_samples samples
    min_samples = 5

    #Multiply all sizes by a scalar
    scaling_factor = 500

    for unique_taxon_name in sorted(unique_taxa_names):
        labels.append(unique_taxon_name)

        for compartment in ['all','mucus','tissue','skeleton']:
            if unique_taxon_name not in compartment_results_dict[compartment].keys():
                continue
            slope = float(compartment_results_dict[compartment][unique_taxon_name]['slope'])
            positive_slope = slope > 0


            result = compartment_results_dict[compartment][unique_taxon_name]['R2']

            if not result or np.isnan(result) or str(result) == 'nan':
                continue

            n_samples = compartment_results_dict[compartment][unique_taxon_name]['N_unique_samples']

            if int(n_samples) < min_samples:
                print("Skipping_taxon:",unique_taxon_name)
                continue


            pval = compartment_results_dict[compartment][unique_taxon_name]['p']
            x_value = x_position_mapping[unique_taxon_name]
            x_values.append(x_value)
            y_values.append(y_position_mapping[compartment]*ytick_scalar)
            sizes.append(result*scaling_factor)
            sig.append(compartment_results_dict[compartment][unique_taxon_name]['p'])
            if positive_slope:
                marker = "^"
            else:
                marker = "v"
            markers.append(marker)

            if float(pval)  < 0.05:
                colors.append(color_mapping[compartment])
                #print(f"{unique_taxon_name} in {compartment}")
                #print("N samples:",n_samples)
                #print('Pval:',pval)
                #print('R2:',result)

                #print("Slope:",slope)
                #print("positive slope?", positive_slope)
                #print("marker",marker)

            else:
                colors.append('whitesmoke')
    return x_values,y_values,colors,sizes,markers,labels,\
      compartments,x_position_mapping,y_position_mapping,color_mapping





if __name__ == "__main__":
    main()
