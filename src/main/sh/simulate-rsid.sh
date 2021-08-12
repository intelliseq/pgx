RSID=$1


SIMULATION_BED="${PROJECT_DIR}resources/miscellaneous/simulated-samples/bed_files/cyps_2mln.bed"
INPUT_TEMPLATE="${PROJECT_DIR}resources/miscellaneous/simulated-samples/samples/pgx-test/create-inputs/template-inputs.json"
CHROMOSOME="chr22"
SAMPLE_NAME="$GENE-star-$STAR"
SIMULATION_VCF_GZ="$OUTPUT_PATH/$GENE-star-$STAR-and-cyps.vcf.gz"
COVERAGE="30"


sed -e "s#\$SIMULATION_BED#$SIMULATION_BED#" \
-e "s#\$CHROMOSOME#$CHROMOSOME#" \
-e "s#\$SAMPLE_NAME#$SAMPLE_NAME#" \
-e "s#\$SIMULATION_VCF_GZ#$SIMULATION_VCF_GZ#" \
-e "s#\$COVERAGE#$COVERAGE#" \
$INPUT_TEMPLATE > $OUTPUT_PATH/$GENE-star-$STAR-simulate-neat-input.json

echo "$OUTPUT_PATH/$GENE-star-$STAR-simulate-neat-input.json"
cat "$OUTPUT_PATH/$GENE-star-$STAR-simulate-neat-input.json"

