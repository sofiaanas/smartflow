SAMPLES = ["AFR","AMR","EAS","FIN","NFE","SAS"]
REPETITION_NUMBER = 1
REPETITION = list(range(1,REPETITION_NUMBER+1))

rule all:
    input:
        expand("sample_pop/{sample}{repetition}.vcf.gz", sample=SAMPLES, repetition=REPETITION),
        "sample_pop/pop_merge.vcf.gz",
        "sample_pop/pop_merge_newid.vcf.gz",
        "sample_pop/pop_merge_newid_corref.vcf.gz",
        "plink/ld_prune.prune.in",
        "plink/ld_prune.prune.out",
        "plink/ld_prune.nosex",
        "plink/ld_prune.log",
        "sample_pop/pop_merge_newid_corref_ldpruned.vcf",
        "smartpca/genotype_rmheader_conum.eigenstratgeno",
        "smartpca/pop.snp",
        "smartpca/pop.ind"


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
        exacindex="databases/ExAC.r0.3.1.sites.vep.reheader.vcf.gz.tbi",
        jar= "../simdrom/simdrom-cli/target/simdrom-cli-0.0.3-SNAPSHOT.jar"
    output:
        "sample_pop/{sample}{repetition}.vcf.gz"
    params:
        sample="{sample}",
        rep="{repetition}"
    shell:
        "java -jar {input.jar} -b {input.exac} -bAC AC_{params.sample} -bAN AN_{params.sample} -n {params.sample}{params.rep} --output {output}"

# Merge all samples into one file.
rule merge:
    input:
        expand("sample_pop/{sample}{repetition}.vcf.gz", sample=SAMPLES, repetition=REPETITION)
    output:
        "sample_pop/pop_merge.vcf.gz"
    shell:
        "bcftools merge {input} -O z -o {output}"

# Give the genes some new IDs (some don't have an ID from ExAC),they will be needed later.
rule new_id:
    input:
        "sample_pop/pop_merge.vcf.gz"
    output:
        "sample_pop/pop_merge_newid.vcf.gz"
    shell:
        """
        zcat {input} | awk -v OFS='\t' '{{ if (/^#/) {{ print $0; next }} else {{ $3=NR; print $0 }} }}' | bgzip -c > {output}
        """

# Correct genotype references.
rule corr_ref:
    input:
        "sample_pop/pop_merge_newid.vcf.gz"
    output:
        "sample_pop/pop_merge_newid_corref.vcf.gz"
    shell:
        "bcftools plugin missing2ref {input} | bgzip -c > {output}"


# Prune regions that display linkage disequilibrium using PLINK.
rule ld_prune:
    input:
        "sample_pop/pop_merge_newid_corref.vcf.gz"
    output:
        "plink/ld_prune.prune.in", # several diffrent filetypes...?
        "plink/ld_prune.prune.out",
        "plink/ld_prune.nosex",
        "plink/ld_prune.log"
    shell:
        "~/PLINK/plink --vcf {input} --indep-pairwise 50 5 0.8 --out plink/ld_prune"

rule rm_regions:
    input:
        prune="plink/ld_prune.prune.out",
        vcf="sample_pop/pop_merge_newid_corref.vcf.gz"
    output:
        uncompressed=temp("sample_pop/pop_merge_newid_corref.vcf"),
        vcf="sample_pop/pop_merge_newid_corref_ldpruned.vcf"
    shell:
        """
        bgzip -dc {input.vcf} > {output.uncompressed};
        awk 'FNR==NR {{a[$i]; next}}; !($3 in a)' {input.prune} {output.uncompressed} > {output}
        """


# Prepare files for smartpca.

# Genotype file: Save only genotype columns, remove header and convert genotypes to suit smartpca.
rule genotype_file:
    input:
        "sample_pop/pop_merge_newid_corref_ldpruned.vcf"
    output:
        "smartpca/genotype_rmheader_conum.eigenstratgeno"
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
                if genotype == "0/0":
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
rule snp_file:
    input:
        "sample_pop/pop_merge_newid_corref_ldpruned.vcf"
    output:
        "smartpca/pop.snp"
    shell:
        """
        cat {input} | grep -v "#" | awk -v OFS='\\t' '{{print $3,$1,"0.0",$2,$4,$5}}' > {output}
        """

# Ind file
rule ind_file:
    input:
    output:
        "smartpca/pop.ind"
    run:
        f = open(output[0],'w')
        for pop in SAMPLES:
            for i in REPETITION:
                print("%s\tU\t%s" % (pop+"sample"+str(i),pop) ,file=f)
        f.close

# Run smartpca
rule smartpca:
    input:
        "smartpca/parfile"
    output:
    shell:
