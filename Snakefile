"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary,
            save_pinned_env

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        env='environment_pinned.yml',
        get_mut_bind_mono_expr=config['mut_bind-monomer_expr'],
        get_mut_bind_di_expr=config['mut_bind-dimer_expr'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        fit_deer_ACE2_titrations='results/summary/compute_deer-ACE2_Kd.md',
        fit_hamster_ACE2_titrations='results/summary/compute_hamster-ACE2_Kd.md',
        fit_bat_ACE2_titrations='results/summary/compute_bat-ACE2_Kd.md',
        fit_cat_ACE2_titrations='results/summary/compute_cat-ACE2_Kd.md',
        deer_ACE2_Kds_file=config['deer_Kds_file'],
        hamster_ACE2_Kds_file=config['hamster_Kds_file'],
        bat_ACE2_Kds_file=config['bat_Kds_file'],
        cat_ACE2_Kds_file=config['cat_Kds_file'],
        collapse_scores='results/summary/collapse_scores.md',
        mut_phenos_file=config['final_variant_scores_mut_file'],
        heatmap_viz=os.path.join(config['visualization_dir'], "heatmap.html")
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:
            
            1. Download prior data for barcode-variant lookup tables and huACE2 binding from prior repositories [1](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS) and [2](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron).
            
            2. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            3. [Fit deer ACE2 titration curves]({path(input.fit_deer_ACE2_titrations)}) to calculate per-barcode K<sub>D,app</sub>, recorded in [this file]({path(input.deer_APN_Kds_file)}).
            
            4. [Fit hamster ACE2 titration curves]({path(input.fit_hamster_ACE2_titrations)}) to calculate per-barcode K<sub>D,app</sub>, recorded in [this file]({path(input.hamster_APN_Kds_file)}).
            
            5. [Fit bat ACE2 titration curves]({path(input.fit_bat_ACE2_titrations)}) to calculate per-barcode K<sub>D,app</sub>, recorded in [this file]({path(input.bat_APN_Kds_file)}).
            
            6. [Fit cat ACE2 titration curves]({path(input.fit_cat_ACE2_titrations)}) to calculate per-barcode K<sub>D,app</sub>, recorded in [this file]({path(input.cat_APN_Kds_file)}).
            
            7. [Derive final genotype-level phenotypes from replicate barcoded sequences]({path(input.collapse_scores)}).
               Generates final phenotypes, recorded in [this file]({path(input.mut_phenos_file)}).

            8. Make interactive data visualizations, available [here](https://tstarrlab.github.io/SARS-CoV-2-RBD_DMS_animal-ACE2/)

            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"


rule save_pinned_env:
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
    log:
    	"environment_pinned.yml"
    shell:
        """
        conda env export > {log}
        """

rule interactive_heatmap_plot:
    """ Make the interactive heatmap for expression and binding.
    """
    input: 
        scores=config['final_variant_scores_mut_file']
    params:
        annotations=config['RBD_sites']
    output:
        html=os.path.join(config['visualization_dir'], "heatmap.html")
    notebook: "RBD-Heatmaps-Interactive-Visualization.ipynb"



rule collapse_scores:
    input:
        config['deer_Kds_file'],
        config['hamster_Kds_file'],
        config['bat_Kds_file'],
        config['cat_Kds_file'],
        config['mut_bind-dimer_expr'],
        config['mut_bind-monomer_expr']
    output:
        config['final_variant_scores_mut_file'],
        md='results/summary/collapse_scores.md',
        md_files=directory('results/summary/collapse_scores_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='collapse_scores.Rmd',
        md='collapse_scores.md',
        md_files='collapse_scores_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_deer_ACE2_titrations:
    input:
        config['codon_variant_table'],
        config['variant_counts_file']
    output:
        config['deer_Kds_file'],
        md='results/summary/compute_deer-ACE2_Kd.md',
        md_files=directory('results/summary/compute_deer-ACE2_Kd_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_deer-ACE2_Kd.Rmd',
        md='compute_deer-ACE2_Kd.md',
        md_files='compute_deer-ACE2_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_hamster_ACE2_titrations:
    input:
        config['codon_variant_table'],
        config['variant_counts_file']
    output:
        config['hamster_Kds_file'],
        md='results/summary/compute_hamster-ACE2_Kd.md',
        md_files=directory('results/summary/compute_hamster-ACE2_Kd_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_hamster-ACE2_Kd.Rmd',
        md='compute_hamster-ACE2_Kd.md',
        md_files='compute_hamster-ACE2_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_bat_ACE2_titrations:
    input:
        config['codon_variant_table'],
        config['variant_counts_file']
    output:
        config['bat_Kds_file'],
        md='results/summary/compute_bat-ACE2_Kd.md',
        md_files=directory('results/summary/compute_bat-ACE2_Kd_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_bat-ACE2_Kd.Rmd',
        md='compute_bat-ACE2_Kd.md',
        md_files='compute_bat-ACE2_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule fit_cat_ACE2_titrations:
    input:
        config['codon_variant_table'],
        config['variant_counts_file']
    output:
        config['cat_Kds_file'],
        md='results/summary/compute_cat-ACE2_Kd.md',
        md_files=directory('results/summary/compute_cat-ACE2_Kd_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_cat-ACE2_Kd.Rmd',
        md='compute_cat-ACE2_Kd.md',
        md_files='compute_cat-ACE2_Kd_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['codon_variant_table'],
        config['barcode_runs']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule get_mut_bind_monomer_expr:
    """Download SARS-CoV-2 DMS for monomer ACE2-binding and expression from URL."""
    output:
        file=config['mut_bind-monomer_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind-monomer_expr_url'], output.file)
        
rule get_mut_bind_dimer_expr:
    """Download SARS-CoV-2 DMS for dimer ACE2-binding and expression from URL."""
    output:
        file=config['mut_bind-dimer_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind-dimer_expr_url'], output.file)
        