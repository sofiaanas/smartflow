SAMPLES = ["AFR","AMR","EAS","FIN","NFE","SAS"]
REPETITION_NUMBER = 1
REPETITION = list(range(1,REPETITION_NUMBER+1))

rule all:
    input:
        expand("sample_pop/{sample}.rep{repetition}.vcf.gz", sample=SAMPLES, repetition=REPETITION),
        "sample_pop/pop_merge.vcf.gz",
        "sample_pop/pop_merge_newid.vcf.gz",
        "sample_pop/pop_merge_newid_corref.vcf.gz",
        "plink/ld_prune",
        "sample_pop/pop_merge_newid_corref_ldpruned.vcf",
        "smartpca/genotype.eigenstratgeno",
        "smartpca/genotype_rmheader.eigenstratgeno",
        "smartpca/genotype_rmheader_conum.eigenstratgeno"


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

rule index_reheader:
    input:
        "databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz"
    output:
        "databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz.tbi"
    shell:
        "tabix {input}"

# Sample individuals of different ethnicities using SIMdrom.
rule sample_pop:
    input:
        exac="databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz",
        exacindex="databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz.tbi"
    output:
        "sample_pop/{sample}.rep{repetition}.vcf.gz"
    params:
        sample="{sample}"
    shell:
        "java -jar simdrom-cli-0.0.1.jar -b {input.exac} -bAC AC_{params.sample} -bAN AN_{params.sample} -n {params.sample} --output {output}"

# Merge all samples into one.
rule merge:
    input:
        expand("sample_pop/{sample}.rep{repetition}.vcf.gz", sample=SAMPLES, repetition=REPETITION)
    output:
        "sample_pop/pop_merge.vcf.gz"
    shell:
        "bcftools merge {input} -O b -o {output}"


# Give the gens new IDs (some don't have an ID from ExAC),they will be needed later.
rule new_id:
    input:
        "sample_pop/pop_merge.vcf.gz"
    output:
        "sample_pop/pop_merge_newid.vcf.gz"
    shell:
        "zcat {input} | awk -v OFS='\t' '{if (NR>226) $3=(NR-226); print $0}'| bgzip -c > {output}"


# Correct genotype references.
rule corr_ref:
    input:
        "sample_pop/pop_merge_newid.vcf.gz"
    output:
        "sample_pop/pop_merge_newid_corref.vcf.gz"
    shell:
        "bcftools plugin missing2ref {input}"


# Prune regions that display linkage disequilibrium using PLINK.
rule ld_prune:
    input:
        "sample_pop/pop_merge_newid_corref.vcf.gz"
    output:
        "plink/ld_prune" # several diffrent filetypes...?
    shell:
        "plink --vcf {input} --indep-pairwise 50 5 0.8 --out {output}"

rule rm_regions:
    input:
        prune="plink/ld_prune.prune.out",
        vcf="sample_pop/pop_merge_newid_corref.vcf.gz"
    output:
        "sample_pop/pop_merge_newid_corref_ldpruned.vcf"
    shell:
        "zcat {input.vcf} | awk 'FNR==NR {a[$i]; next}; !($3 in a)' {input.prune} > {output}"


# Prepare files for smartpca.

# Genotype file:
# Save only genotype columns.
rule gen_col:
    input:
        "sample_pop/pop_merge_newid_corref_ldpruned.vcf"
    output:
        "smartpca/genotype.eigenstratgeno"
    shell:
        "cut -f10-27 {input} > {output}"

# Remove header.
rule rm_header:
    input:
        "smartpca/genotype.eigenstratgeno"
    output:
        "smartpca/genotype_rmheader.eigenstratgeno"
    shell:
        "awk 'NR > 226' {input} > {output}"

# Convert genotypes to suit smartpca.
rule con_num_2:
    input:
        "smartpca/genotype_rmheader.eigenstratgeno"
    output:
        "smartpca/genotype_rmheader_conum.eigenstratgeno"
    shell:
        "sed -i 's:0/0:2:g' {input} > {output}"
#awk '{if ($0==0/0) $0=2; print $0}'
#sed -i 's:0/0:2:g'



# SNP file
# Ind file

# Run smartpca
