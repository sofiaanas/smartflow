
CHROMOSOMES=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT"]
CHROMEX=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
AUTOSOMES=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
SAMPLES = ["AFR","AMR","EAS","FIN","NFE","SAS"]
REAL = ["NA19020","HG01550","NA18620","HG00310","HG00140","HG02600"]

REPETITION_NUMBER = 10
REPETITION = list(range(1,REPETITION_NUMBER+1))

ANALYSIS = ["merge_exac_af","merge_exac_ac"]

rule all:
    input:
        expand("region/Panel.{chex}.coverage.txt.gz", chex=CHROMEX),
        expand("region/Panel.{chex}.coverage.txt.gz.tbi", chex=CHROMEX),
        "region/regions.bed",
        expand("1KG/samples/1KGsamples.{auto}.vcf.gz", auto=AUTOSOMES),
        expand("1KG/samples/1KGsamples.{auto}.vcf.gz.tbi", auto=AUTOSOMES),
        "1KG/samples/1KGsamples.chrMT.vcf.gz",
        "1KG/samples/1KGsamples.chrX.vcf.gz",
        "1KG/samples/1KGsamples.chrY.vcf.gz",
        "1KG/samples/1KGsamples.chrMT.vcf.gz.tbi",
        "1KG/samples/1KGsamples.chrX.vcf.gz.tbi",
        "1KG/samples/1KGsamples.chrY.vcf.gz.tbi",
        "1KG/1KGsamples_concat.vcf.gz",
        "1KG/1KGsamples_concat.vcf.gz.tbi",
        "1KG/1KG.regcut.vcf.gz",
        "1KG/1KG.regcut.vcf.gz.tbi",
        "ExAC/ExAC.regcut.vcf.gz",
        "ExAC/ExAC.regcut.vcf.gz.tbi",
        # samplings
        expand("samplings/ExAC/sample_AF_pop/{sample}.{repetition}.vcf.gz", sample=SAMPLES, repetition=REPETITION),
        # analysis
        expand("analysis/{analysis}/prepsmart/pop_merge.vcf.gz", analysis=ANALYSIS),
        expand("analysis/{analysis}/prepsmart/pop_merge_newid.vcf.gz", analysis=ANALYSIS),
        expand("analysis/{analysis}/prepsmart/pop_merge_newid_corref.vcf.gz", analysis=ANALYSIS),
        expand("analysis/{analysis}/plink/ld_prune.prune.in", analysis=ANALYSIS),
        expand("analysis/{analysis}/plink/ld_prune.prune.out", analysis=ANALYSIS),
        expand("analysis/{analysis}/plink/ld_prune.nosex", analysis=ANALYSIS),
        expand("analysis/{analysis}/plink/ld_prune.log", analysis=ANALYSIS),
        expand("analysis/{analysis}/prepsmart/pop_merge_newid_corref_ldpruned.vcf", analysis=ANALYSIS),
        expand("analysis/{analysis}/smartpca/genfile.eigenstratgeno", analysis=ANALYSIS),
        expand("analysis/{analysis}/smartpca/snpfile.snp", analysis=ANALYSIS),
        expand("analysis/{analysis}/smartpca/indfile.ind", analysis=ANALYSIS),
        expand("analysis/{analysis}/smartpca/pca.log", analysis=ANALYSIS),
        expand("analysis/{analysis}/smartpca/pca.evec", analysis=ANALYSIS),
        expand("analysis/{analysis}/smartpca/pca.eval", analysis=ANALYSIS),
        # "smartpca/grmjunk",
        # "smartpca/grmjunk.id",
        expand("analysis/{analysis}/smartpca/plotpca.xtxt", analysis=ANALYSIS),
        expand("analysis/{analysis}/smartpca/plotpca.ps", analysis=ANALYSIS)
        # "smartpca/plotpca.pdf"


# Download the regions for the exome and their indexes
rule download_region:
    input:
    output:
        expand("region/Panel.{chex}.coverage.txt.gz", chex=CHROMEX)
    run:
        for chr in CHROMEX:
            shell("wget -O region/Panel.{chr}.coverage.txt.gz ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/coverage/Panel.{chr}.coverage.txt.gz")

rule download_regionidx:
    input:
    output:
        expand("region/Panel.{chex}.coverage.txt.gz.tbi", chex=CHROMEX)
    run:
        for chr in CHROMEX:
            shell("wget -O region/Panel.{chr}.coverage.txt.gz.tbi ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/coverage/Panel.{chr}.coverage.txt.gz.tbi")

# Pick regions with given covrage
rule pick_region:
    input:
        dir="region/",
        scr="scripts/regions.py"
    output:
        "region/regions.bed"
    shell:
        "python {input.scr} --coverage 20 --directory {input.dir} > {output}"


####################################  1KG  ####################################


# Download the real genomes and indexes from 1000 Genomes project, first the autosomes
rule download_auto:
    input:
    output:
        temp=(expand("1KG/download/ALL.{auto}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", auto=AUTOSOMES))
    run:
        for chr in AUTOSOMES:
            shell("wget -O 1KG/download/ALL.{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")


rule download_autoidx:
    input:
    output:
        temp=(expand("1KG/download/ALL.{auto}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi", auto=AUTOSOMES))
    run:
        for chr in AUTOSOMES:
            shell("wget -O 1KG/download/ALL.{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi")

# Pick out six samples from each autosome (one from each population) defined in samples.txt
rule pick_auto:
    input:
        samples="samples.txt",
        vcf="1KG/download/ALL.{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        index="1KG/download/ALL.{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi"
    output:
        "1KG/samples/1KGsamples.{chr}.vcf.gz"
    shell:
        "bcftools view -Oz --min-ac=1 -S {input.samples} {input.vcf} > {output}"

# Create new indexes
rule autosample_idx:
	input:
		"1KG/samples/1KGsamples.{chr}.vcf.gz"
	output:
		"1KG/samples/1KGsamples.{chr}.vcf.gz.tbi"
	shell:
		"tabix {input}"

# Download the real genomes and indexes from 1000 Genomes project, now the gonosomes
rule download_gono:
	input:
	output:
		mt=temp("1KG/download/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"),
		x=temp("1KG/download/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"),
		y=temp("1KG/download/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz")
	shell:
		"wget -O {output.mt} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz && "
		"wget -O {output.x} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz && "
		"wget -O {output.y} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz "

rule download_gonoidx:
	input:
	output:
		mt=temp("1KG/download/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi"),
		x=temp("1KG/download/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi"),
		y=temp("1KG/download/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi")
	shell:
		"wget -O {output.mt} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi && "
		"wget -O {output.x} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi && "
		"wget -O {output.y} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi "

# Pick out six samples from each gonosome (one from each population) defined in samples.txt
rule pick_gono:
    input:
        samples="samples.txt",
        mt="1KG/download/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz",
        mtidx="1KG/download/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi",
        x="1KG/download/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz",
        xidx="1KG/download/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi",
        y="1KG/download/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz",
        yidx="1KG/download/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz.tbi"
    output:
        mt="1KG/samples/1KGsamples.chrMT.vcf.gz",
        x="1KG/samples/1KGsamples.chrX.vcf.gz",
        y="1KG/samples/1KGsamples.chrY.vcf.gz"
    shell:
        """
        bcftools view -Oz --min-ac=1 -S {input.samples} {input.mt} > {output.mt};
        bcftools view -Oz --min-ac=1 -S {input.samples} {input.x} > {output.x};
        bcftools view -Oz --min-ac=1 -S {input.samples} {input.y} > {output.y}
        """

# Create new indexes
rule gonosample_idx:
    input:
        mt="1KG/samples/1KGsamples.chrMT.vcf.gz",
        x="1KG/samples/1KGsamples.chrX.vcf.gz",
        y="1KG/samples/1KGsamples.chrY.vcf.gz"
    output:
        mt="1KG/samples/1KGsamples.chrMT.vcf.gz.tbi",
        x="1KG/samples/1KGsamples.chrX.vcf.gz.tbi",
        y="1KG/samples/1KGsamples.chrY.vcf.gz.tbi"
    shell:
        """
        tabix {input.mt};
        tabix {input.x};
        tabix {input.y}
        """

# Concatenate the 1KG chromosome vcfs into one vcf
rule concat1KG:
	input:
		vcf=expand("1KG/samples/1KGsamples.{chrom}.vcf.gz", chrom=CHROMOSOMES),
		index=expand("1KG/samples/1KGsamples.{chrom}.vcf.gz.tbi", chrom=CHROMOSOMES)
	output:
		"1KG/1KGsamples_concat.vcf.gz"
	shell:
		"bcftools concat {input.vcf} | bcftools view -v snps,indels | bgzip -c > {output}"

# Create index
rule concat_idx:
	input:
		"1KG/1KGsamples_concat.vcf.gz"
	output:
		"1KG/1KGsamples_concat.vcf.gz.tbi"
	shell:
		"tabix {input}"

rule filter_1kg:
    input:
        vcf="1KG/1KGsamples_concat.vcf.gz",
        bed="region/regions.bed"
    output:
        "1KG/1KG.regcut.vcf.gz"
    shell:
        "bedtools intersect -header -a {input.vcf} -b {input.bed} 2>/dev/null | bgzip -c > {output}"

rule filterkg_idx:
    input:
        "1KG/1KG.regcut.vcf.gz"
    output:
        "1KG/1KG.regcut.vcf.gz.tbi"
    shell:
        "tabix {input}"

###################################  EXAC  ###################################


# Download and prepare ExAC files for SIMdrom.
rule download_exac:
    input:
    output:
        vcf="databases/ExAC.r0.3.1.sites.vep.vcf.gz",
        index="databases/ExAC.r0.3.1.sites.vep.vcf.gz.tbi"
    shell:
        """
        wget -O {output.vcf} ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz;
        wget -O {output.index} ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz.tbi
        """

rule reheader_exac:
    input:
        corr_h="files/correct_exac_header.txt",
        exac="databases/ExAC.r0.3.1.sites.vep.vcf.gz"
    output:
        "databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz"
    shell:
        "tabix -r {input.corr_h} {input.exac} > {output}"

rule reheader_idx:
    input:
        "databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz"
    output:
        "databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz.tbi"
    shell:
        "tabix {input}"

rule filter_exac:
    input:
        vcf="databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz",
        bed="region/regions.bed"
    output:
        "ExAC/ExAC.regcut.vcf.gz"
    shell:
        "bedtools intersect -header -a {input.vcf} -b {input.bed} 2>/dev/null | bgzip -c > {output}"

rule filterex_idx:
    input:
        "ExAC/ExAC.regcut.vcf.gz"
    output:
        "ExAC/ExAC.regcut.vcf.gz.tbi"
    shell:
        "tabix {input}"


############################################################ Sampling

# Sample individuals of different ethnicities using SIMdrom.
rule sample_AF_pop:
    input:
        exac="ExAC/ExAC.regcut.vcf.gz",
        eidx="ExAC/ExAC.regcut.vcf.gz.tbi",
        jar= "../simdrom/simdrom-cli/target/simdrom-cli-0.0.3-SNAPSHOT.jar"
    output:
        "samplings/ExAC/sample_AF_pop/{sample}.{repetition}.vcf.gz"
    params:
        sample="{sample}",
        rep="{repetition}"
    shell:
        "java -jar {input.jar} -b {input.exac} -bAC AC_{params.sample} -bAN AN_{params.sample} -n {params.sample}.{params.rep} --output {output}"

rule sample_AC_pop:
    input:
        exac="ExAC/ExAC.regcut.vcf.gz",
        eidx="ExAC/ExAC.regcut.vcf.gz.tbi",
        jar= "../simdrom/simdrom-cli/target/simdrom-cli-0.0.3-SNAPSHOT.jar"
    output:
        "samplings/ExAC/sample_AC_pop/{sample}.{repetition}.vcf.gz"
    params:
        sample="{sample}",
        rep="{repetition}"
    shell:
        "java -jar {input.jar} -b {input.exac} -bAChom Hom_{params.sample} -bAChet Het_{params.sample} -bAChemi Hemi_{params.sample} -bAN AN_{params.sample} -n {params.sample}.{params.rep} --output {output}"

############################ analysis

# Merge all samples (exac and real) into one file.
rule merge_exac_af:
    input:
        exac=expand("samplings/ExAC/sample_AF_pop/{sample}.{repetition}.vcf.gz", sample=SAMPLES, repetition=REPETITION),
        kg="1KG/1KG.regcut.vcf.gz"
    output:
        "analysis/merge_exac_af/prepsmart/pop_merge.vcf.gz"
    shell:
        "bcftools merge {input.exac} {input.kg} -O z -o {output}"

rule merge_exac_ac:
    input:
        exac=expand("samplings/ExAC/sample_AC_pop/{sample}.{repetition}.vcf.gz", sample=SAMPLES, repetition=REPETITION),
        kg="1KG/1KG.regcut.vcf.gz"
    output:
        "analysis/merge_exac_ac/prepsmart/pop_merge.vcf.gz"
    shell:
        "bcftools merge {input.exac} {input.kg} -O z -o {output}"

###############################  prep smartpca  ###############################


# Give the snps new IDs (some don't have an ID from ExAC),they will be needed later.
rule new_id:
    input:
        "analysis/{analysis}/prepsmart/pop_merge.vcf.gz"
    output:
        "analysis/{analysis}/prepsmart/pop_merge_newid.vcf.gz"
    shell:
        """
        zcat {input} | awk -v OFS='\t' '{{ if (/^#/) {{ print $0; next }} else {{ $3=NR; print $0 }} }}' | bgzip -c > {output}
        """

# Correct the genotype references.
rule corr_ref:
    input:
        "analysis/{analysis}/prepsmart/pop_merge_newid.vcf.gz"
    output:
        "analysis/{analysis}/prepsmart/pop_merge_newid_corref.vcf.gz"
    shell:
        "bcftools plugin missing2ref {input} | bgzip -c > {output}"

# Prune regions that display linkage disequilibrium using PLINK.
rule ld_prune:
    input:
        "analysis/{analysis}/prepsmart/pop_merge_newid_corref.vcf.gz"
    output:
        "analysis/{analysis}/plink/ld_prune.prune.in",
        "analysis/{analysis}/plink/ld_prune.prune.out",
        "analysis/{analysis}/plink/ld_prune.nosex",
        "analysis/{analysis}/plink/ld_prune.log"
    params:
        analysis="{analysis}"
    shell:
        "~/PLINK/plink --vcf {input} --indep-pairwise 50 5 0.8 --out analysis/{params.analysis}/plink/ld_prune"

rule rm_regions:
    input:
        prune="analysis/{analysis}/plink/ld_prune.prune.out",
        vcf="analysis/{analysis}/prepsmart/pop_merge_newid_corref.vcf.gz"
    output:
        uncompressed=temp("analysis/{analysis}/prepsmart/pop_merge_newid_corref.vcf"),
        vcf="analysis/{analysis}/prepsmart/pop_merge_newid_corref_ldpruned.vcf"
    shell:
        """
        bgzip -dc {input.vcf} > {output.uncompressed};
        awk 'FNR==NR {{a[$i]; next}}; !($3 in a)' {input.prune} {output.uncompressed} > {output.vcf};
        """

# Genotype file: Save only genotype columns, remove header and convert genotypes to suit smartpca.
rule genotypefile:
    input:
        "analysis/{analysis}/prepsmart/pop_merge_newid_corref_ldpruned.vcf"
    output:
        "analysis/{analysis}/smartpca/genfile.eigenstratgeno"
    run:
        fin = open(input[0], 'r')
        fout = open(output[0],'w')
        for line in fin:
            if line.startswith("#"):
                continue
            split_line = line.split("\t")
            genotypes = split_line[9:len(split_line)]
            results = []
            for genotype in genotypes:
                genotype = genotype.strip()
                if (genotype == "0/0" or genotype == "0|0"):
                    result = "2"
                elif "0" in genotype:
                    result = "1"
                else:
                    result = "0"
                results.append(result)
            print("".join(results), file=fout)

        fin.close()
        fout.close()

# SNP file
rule snpfile_exac:
    input:
        "analysis/{analysis}/prepsmart/pop_merge_newid_corref_ldpruned.vcf"
    output:
        "analysis/{analysis}/smartpca/snpfile.snp"
    shell:
        """
        cat {input} | grep -v "#" | awk -v OFS='\\t' '{{ if ($1 == "X") $1=23; else if ($1=="Y") $1=24; print $3,$1,"0.0",$2,$4,$5}}' > {output}
        """

# Individuals file
rule indfile:
    input:
        "analysis/{analysis}/prepsmart/pop_merge_newid_corref_ldpruned.vcf"
    output:
        "analysis/{analysis}/smartpca/indfile.ind"
    shell:
        """
        bcftools query -l {input} | awk -v OFS="\t" -F"." '{{print $1$2,"U",$1}}' > {output}
        """

# Create the parfile
rule smartpcaconfig:
    input:
        genotypename = "analysis/{analysis}/smartpca/genfile.eigenstratgeno",
        snpname="analysis/{analysis}/smartpca/snpfile.snp",
        indivname="analysis/{analysis}/smartpca/indfile.ind"
    output:
        "analysis/{analysis}/smartpca/parfile"
    params:
        evecoutname="analysis/{analysis}/smartpca/pca.evec",
        evaloutname="analysis/{analysis}/smartpca/pca.eval"
    shell:
        """
        echo "genotypename: {input.genotypename}" >> {output};
        echo "snpname: {input.snpname}" >> {output};
        echo "indivname: {input.indivname}" >> {output};
        echo "evecoutname: {params.evecoutname}" >> {output};
        echo "evaloutname: {params.evaloutname}" >> {output};
        echo "numoutlieriter:   0" >> {output};
        echo "ealtnormstyle:    NO" >> {output};
        echo "familynames:     NO" >> {output};
        echo "grmoutname:      grmjunk" >> {output};
        """

# Run smartpca
rule smartpca:
    input:
        "analysis/{analysis}/smartpca/parfile"
    output:
        log="analysis/{analysis}/smartpca/pca.log",
        evec="analysis/{analysis}/smartpca/pca.evec",
        evl="analysis/{analysis}/smartpca/pca.eval"
        # grmjunk="smartpca/grmjunk",
        # grmid="smartpca/grmjunk.id"
    shell:
        "../../EIG-6.1.4/bin/smartpca -p {input} > {output.log}"

# Plot the results from smartpca
rule plot:
    input:
        "analysis/{analysis}/smartpca/pca.evec"
    output:
        xtxt="analysis/{analysis}/smartpca/plotpca.xtxt",
        ps="analysis/{analysis}/smartpca/plotpca.ps",
        # pdf="smartpca/plotpca.pdf"

    params:
        sample=":".join(SAMPLES+REAL)
    shell:
        "../../EIG-6.1.4/bin/ploteig -i {input} -c 1:2 -p {params.sample} -x -o {output.xtxt}"
