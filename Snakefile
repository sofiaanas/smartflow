SAMPLES = ["AFR","AMR","EAS","FIN","NFE","SAS"]

rule sample_pop:
    input:
        "ExAC.r0.3.1.sites.vep.reheader.vcf.gz"
    output:
        "sample_pop/{sample}.vcf.gz"
    shell:
        "java -jar simdrom-cli-0.0.1.jar -b {input} -bAC AC_{sample} -bAN AN_{sample} -n {sample} --output {output}"
