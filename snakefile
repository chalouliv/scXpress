
configfile:"./config.yaml"
# -----------------------------------------------
# Input function to generate rule all files 
# -----------------------------------------------

def generate_output_files(config):
    outs = []
    for ref in config['signatures']:
       if (config['signatures'][ref]['path']=="msigdbr"):
            temp = expand("{out}/{query}/{signature}/{query}.ex.{signature}.{method}.csv",
                    out=config["out"],
                    query=config["queries"],
                    signature=ref,
                    method=config["methods"])
       else:
            temp = expand("{out}/{query}/{signature}/{query}.in.{signature}.{method}.csv",
                    out=config["out"],
                    query=config["queries"],
                    signature=ref,
                    method=config["methods"])
       outs += temp
    return outs


#def generate_q_paths(config):
#    q_paths = ""
#    for query in config['queries']:
#        temp = config['queries'][query]['path'] + " "
#        print(temp)
#	q_paths += temp
#    return q_paths

#def generate_m_paths(config):
#    m_paths = ""
#    for query in config['queries']:
#        temp = config['queries'][query]['meta'] + " "
#        print(temp)
#	m_paths += temp
#    return m_paths

# -----------------------------------------------
# rule all
# -----------------------------------------------

rule all:
    input: 
    	file = expand("{out}/{report}",
    	  out = config["out"],
        report = config["report"]+ ".html")
       
# -----------------------------------------------
# GSVA 
# -----------------------------------------------

rule GSVA_in:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
        signature_path = lambda wildcards: config['signatures'][wildcards.signature]['path']	
    log:
        "{out}/logs/{query}.in.{signature}.gsva.log"
    output:
        "{out}/{query}/{signature}/{query}.in.{signature}.gsva.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        specifier = lambda wildcards: config['signatures'][wildcards.signature]['specifier']
    shell:
        """
        Rscript scripts/gsva_run.R \
            {output} \
            {input.query_path} \
            {input.signature_path} \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.specifier} \
            gsva \
	          &> {log}
        """
        
rule GSVA_ex:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
    log:
        "{out}/logs/{query}.ex.{signature}.gsva.log"
    output:
        "{out}/{query}/{signature}/{query}.ex.{signature}.gsva.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        species = lambda wildcards: config['signatures'][wildcards.signature]['species'],
        category = lambda wildcards: config['signatures'][wildcards.signature]['category'],
        subcategory = lambda wildcards: config['signatures'][wildcards.signature]['subcategory'],
        gs_style = lambda wildcards: config['signatures'][wildcards.signature]['gs_style'],
        gene_style = lambda wildcards: config['signatures'][wildcards.signature]['gene_style']
    shell:
        """
        Rscript scripts/gsva_run.R \
            {output} \
            {input.query_path} \
            msigdbr \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.species} \
            {params.category} \
            {params.subcategory} \
            {params.gs_style} \
            {params.gene_style}\
            gsva \
	          &> {log}
        """
        
# -----------------------------------------------
# GSVA w/ Gaussian
# -----------------------------------------------

rule GSVA_gaussian_in:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
        signature_path = lambda wildcards: config['signatures'][wildcards.signature]['path']	
    log:
        "{out}/logs/{query}.in.{signature}.gsva_gaussian.log"
    output:
        "{out}/{query}/{signature}/{query}.in.{signature}.gsva_gaussian.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        specifier = lambda wildcards: config['signatures'][wildcards.signature]['specifier']
    shell:
        """
        Rscript scripts/gsva_run.R \
            {output} \
            {input.query_path} \
            {input.signature_path} \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.specifier} \
            gsva_gaussian \
	          &> {log}
        """
        
rule GSVA_gaussian_ex:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
    log:
        "{out}/logs/{query}.ex.{signature}.gsva_gaussian.log"
    output:
        "{out}/{query}/{signature}/{query}.ex.{signature}.gsva_gaussian.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        species = lambda wildcards: config['signatures'][wildcards.signature]['species'],
        category = lambda wildcards: config['signatures'][wildcards.signature]['category'],
        subcategory = lambda wildcards: config['signatures'][wildcards.signature]['subcategory'],
        gs_style = lambda wildcards: config['signatures'][wildcards.signature]['gs_style'],
        gene_style = lambda wildcards: config['signatures'][wildcards.signature]['gene_style']
    shell:
        """
        Rscript scripts/gsva_run.R \
            {output} \
            {input.query_path} \
            msigdbr \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.species} \
            {params.category} \
            {params.subcategory} \
            {params.gs_style} \
            {params.gene_style}\
            gsva_gaussian \
	          &> {log}
        """
        
        
        
# -----------------------------------------------
# SSGSEA
# -----------------------------------------------

rule SSGSEA_in:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
        signature_path = lambda wildcards: config['signatures'][wildcards.signature]['path']	
    log:
        "{out}/logs/{query}.in.{signature}.ssgsea.log"
    output:
        "{out}/{query}/{signature}/{query}.in.{signature}.ssgsea.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        specifier = lambda wildcards: config['signatures'][wildcards.signature]['specifier']
    shell:
        """
        Rscript scripts/gsva_run.R \
            {output} \
            {input.query_path} \
            {input.signature_path} \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.specifier} \
            ssgsea \
	          &> {log}
        """
        
rule SSGSEA_ex:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
    log:
        "{out}/logs/{query}.ex.{signature}.ssgsea.log"
    output:
        "{out}/{query}/{signature}/{query}.ex.{signature}.ssgsea.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        species = lambda wildcards: config['signatures'][wildcards.signature]['species'],
        category = lambda wildcards: config['signatures'][wildcards.signature]['category'],
        subcategory = lambda wildcards: config['signatures'][wildcards.signature]['subcategory'],
        gs_style = lambda wildcards: config['signatures'][wildcards.signature]['gs_style'],
        gene_style = lambda wildcards: config['signatures'][wildcards.signature]['gene_style']
    shell:
        """
        Rscript scripts/gsva_run.R \
            {output} \
            {input.query_path} \
            msigdbr \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.species} \
            {params.category} \
            {params.subcategory} \
            {params.gs_style} \
            {params.gene_style} \
            ssgsea \
            &> {log}
        """

# -----------------------------------------------
# SSGSEA w/ Normalization
# -----------------------------------------------

rule SSGSEA_norm_in:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
        signature_path = lambda wildcards: config['signatures'][wildcards.signature]['path']	
    log:
        "{out}/logs/{query}.in.{signature}.ssgsea_norm.log"
    output:
        "{out}/{query}/{signature}/{query}.in.{signature}.ssgsea_norm.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        specifier = lambda wildcards: config['signatures'][wildcards.signature]['specifier']
    shell:
        """
        Rscript scripts/gsva_run.R \
            {output} \
            {input.query_path} \
            {input.signature_path} \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.specifier} \
            ssgsea_norm \
	          &> {log}
        """
        
rule SSGSEA_norm_ex:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
    log:
        "{out}/logs/{query}.ex.{signature}.ssgsea_norm.log"
    output:
        "{out}/{query}/{signature}/{query}.ex.{signature}.ssgsea_norm.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        species = lambda wildcards: config['signatures'][wildcards.signature]['species'],
        category = lambda wildcards: config['signatures'][wildcards.signature]['category'],
        subcategory = lambda wildcards: config['signatures'][wildcards.signature]['subcategory'],
        gs_style = lambda wildcards: config['signatures'][wildcards.signature]['gs_style'],
        gene_style = lambda wildcards: config['signatures'][wildcards.signature]['gene_style']
    shell:
        """
        Rscript scripts/gsva_run.R \
            {output} \
            {input.query_path} \
            msigdbr \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.species} \
            {params.category} \
            {params.subcategory} \
            {params.gs_style} \
            {params.gene_style} \
            ssgsea_norm \
            &> {log}
        """
        
# -----------------------------------------------
# Module Score 
# -----------------------------------------------

rule modscore_in:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
        signature_path = lambda wildcards: config['signatures'][wildcards.signature]['path']	
    log:
        "{out}/logs/{query}.in.{signature}.module_score.log"
    output:
        "{out}/{query}/{signature}/{query}.in.{signature}.module_score.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        specifier = lambda wildcards: config['signatures'][wildcards.signature]['specifier']
    shell:
        """
        Rscript scripts/modscore_run.R \
            {output} \
            {input.query_path} \
            {input.signature_path}\
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.specifier} \
	          &> {log}
        """
        
rule modscore_ex:
    input:
        query_path = lambda wildcards: config["queries"][wildcards.query]['path'],
    log:
        "{out}/logs/{query}.ex.{signature}.module_score.log"
    output:
        "{out}/{query}/{signature}/{query}.ex.{signature}.module_score.csv"
    params:
        out = config['out'],
        query_meta = lambda wildcards: config["queries"][wildcards.query]['meta'],
        gs_name = lambda wildcards: config['signatures'][wildcards.signature]['gs_name'],
        avg_by = lambda wildcards: config['queries'][wildcards.query]['avg_by'],
        species = lambda wildcards: config['signatures'][wildcards.signature]['species'],
        category = lambda wildcards: config['signatures'][wildcards.signature]['category'],
        subcategory = lambda wildcards: config['signatures'][wildcards.signature]['subcategory'],
        gs_style = lambda wildcards: config['signatures'][wildcards.signature]['gs_style'],
        gene_style = lambda wildcards: config['signatures'][wildcards.signature]['gene_style']
    shell:
        """
        Rscript scripts/modscore_run.R \
            {output} \
            {input.query_path} \
            msigdbr \
            {params.query_meta} \
            {params.avg_by} \
            {params.gs_name} \
            {params.species} \
            {params.category} \
            {params.subcategory} \
            {params.gs_style} \
            {params.gene_style} \
            &> {log}
        """
        
rule Report:
    input:
        file = generate_output_files(config)
    log: "{out}/logs/{report}.log"
    output: "{out}/{report}.html"
    params:
        out = config["out"],
        report = config["report"],
        queries = expand("{query}", query=config["queries"]),
	lists = expand("{list}", list=config["signatures"]),
        paths = expand("{path}", path=[config["signatures"][sig]["path"] for sig in config["signatures"]]),
        methods = expand("{method}", method=config["methods"])
    shell:
        """
        Rscript -e "rmarkdown::render(
            'scripts/my_report.Rmd',
            params = list(
              lists = '{params.lists}',
              paths = '{params.paths}',
              methods = '{params.methods}',
              queries = '{params.queries}',
              out_directory = '{params.out}'
            ),
            output_file = '{params.report}.html',
            output_dir = '{params.out}'
        )" &> {log}
        """
