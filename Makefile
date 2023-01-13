.PHONY: all clean list



FB_REL=6.48
TRANSPOSONS_REL=10.2

REF_GENOME=ref/dmel-all-chromosome-r${FB_REL}.fasta
REF_TRANSPOSONS=ref/dmel-all-transposon-r${FB_REL}.fasta
TRANSPOSONS_DMEL=transposons/D_mel_transposon_sequence_set_v${TRANSPOSONS_REL}.fa


SX4_FOR=1_R1_001 2_R1_001 7_R1_001 8_R1_001
SX4_REV=1_R2_001 2_R2_001 7_R2_001 8_R2_001
SX4=${SX4_FOR} ${SX4_REV}



SX4_COV=${addprefix output/fast_coverage/, ${addsuffix .txt, ${SX4}}}
SX4_TE_MAP=${addprefix output/te_mapper/, ${addsuffix /te_mapper_output.json, ${SX4_FOR} ${SX4_REV}}}
MISC_1=${addprefix output/te_mapper/, ${addsuffix /split_reads_for_tgt.csv, ${SX4_FOR} ${SX4_REV}}}
MISC_2=${addprefix output/te_mapper/, ${addsuffix /figure_SX4Ch7.png, ${SX4_FOR} ${SX4_REV}}}

TARGETS=${SX4_COV} ${SX4_TE_MAP} ${MISC_1} ${MISC_2}

all: ${TARGETS}

clean:
	rm -rf output/*

list:
	@echo ${TARGETS} | xargs -n 1

output/fast_coverage/%.txt:
	@mkdir -p ${dir $@}
	cargo run --bin calculate_stats --release -- fast_coverage -q reads/${basename ${notdir $@}}.fastq -r ${REF_GENOME} -o $@

output/te_mapper/%/te_mapper_output.json: reads/%.fastq
	@mkdir -p ${dir $@}
	sx map -j \
	--reads $< \
	--ref ref/dmel-all-chromosome-r6.48.fasta \
	--transposons transposons/D_mel_transposon_sequence_set_v10.2.fa \
	--result ${dir $@}

# extract target split reads (corresponding to SX-4 insertions into natural TE's,
# confirmed through inverse PCR) into CSV file
output/te_mapper/%/split_reads_for_tgt.csv: scripts/get_table1_split_reads.py output/te_mapper/%/te_mapper_output.json
	python3 $< ${lastword ${subst /, , ${dir $@}}}

# generate TE figures
output/te_mapper/%/figure_SX4Ch7.png: scripts/mk_table1_diagrams.py output/te_mapper/%/split_reads_for_tgt.csv
	python3 $< ${lastword ${subst /, , ${dir $@}}}