SAMPLES = ["AFR","AMR","EAS","FIN","NFE","SAS"]
REPETITION_NUMBER = 1
REPETITION = list(range(1,REPETITION_NUMBER+1))

rule all:
    input:
        expand("sample_pop/{sample}.rep{repetition}.vcf.gz", sample=SAMPLES, repetition=REPETITION)

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
'''
rule remove_regions:
    input:
        "sample_pop/{sample}.rep{repetition}.vcf.gz"
    output:
        "sample_pop/{sample}.rep{repetition}.filtered.vcf.gz"
'''
