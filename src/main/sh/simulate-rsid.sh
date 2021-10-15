RSID=$1
DBSNP="resources/GRCh38.dbSNP155.chr.norm.vcf.gz"
VCF="resources/pgx/vcf/$RSID.vcf.gz"
zcat $DBSNP | grep -P "^#|\t$RSID\t" | bgzip -c > $VCF

SIMULATION_BED="pgx/src/main/resources/bed/cyps_2mln.bed"
INPUT_TEMPLATE="pgx/src/main/resources/json/template-inputs.json"
OPTIONS_TEMPLATE="pgx/src/main/resources/json/template-options.json"
GENO_INPUT_TEMPLATE="pgx/src/main/resources/json/template-genotyping-inputs.json"
CHROMOSOME=$(zcat $VCF | tail -1 | cut -f 1)
SAMPLE_NAME=$RSID
SIMULATION_VCF_GZ=$VCF
COVERAGE="30"
SIMULATED_PATH="resources/pgx/simulated/$RSID"

sed -e "s#\$SIMULATION_BED#$SIMULATION_BED#" \
-e "s#\$CHROMOSOME#$CHROMOSOME#" \
-e "s#\$SAMPLE_NAME#$SAMPLE_NAME#" \
-e "s#\$SIMULATION_VCF_GZ#$SIMULATION_VCF_GZ#" \
-e "s#\$COVERAGE#$COVERAGE#" \
$INPUT_TEMPLATE > resources/pgx/inputs/$RSID-input.json

sed -e "s#\$PATH#$SIMULATED_PATH#" \
$OPTIONS_TEMPLATE > resources/pgx/options/$RSID-options.json

cromwell run pgx/src/main/wdl/simulate-neat.wdl --inputs resources/pgx/inputs/$RSID-input.json --options resources/pgx/options/$RSID-options.json

FQ1=$SIMULATED_PATH"/"$RSID"_1.fq.gz"
FQ2=$SIMULATED_PATH"/"$RSID"_2.fq.gz"
sed -e "s#\$FQ1#$FQ1#" \
    -e "s#\$FQ2#$FQ2#" \
$GENO_INPUT_TEMPLATE > resources/pgx/inputs/$RSID-geno-inputs.json

GENOTYPING_PATH="resources/pgx/genotyping/$RSID"
sed -e "s#\$PATH#$GENOTYPING_PATH#" \
$OPTIONS_TEMPLATE > resources/pgx/options/$RSID-geno-options.json

cromwell run pgx/src/main/wdl/pgx.wdl --inputs resources/pgx/inputs/$RSID-geno-inputs.json --options resources/pgx/options/$RSID-geno-options.json
zcat resources/pgx/genotyping/$RSID/no_id_provided_uniq.g.vcf.gz | grep -P "#|0/1|1/1" | bcftools norm -m- | grep -P "#|0/1|1/1" | grep -vP '\t44057827\t' | bgzip -c > resources/pgx/genotyping/$RSID/filtered.vcf.gz; tabix -p vcf resources/pgx/genotyping/$RSID/filtered.vcf.gz;
bcftools annotate -c ID -a resources/pgx/cyp2d6.rsid.vcf.gz resources/pgx/genotyping/$RSID/filtered.vcf.gz | bgzip -c > resources/pgx/genotyping/$RSID/$RSID-final.vcf.gz
