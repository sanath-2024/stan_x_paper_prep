.PHONY: all clean



FB_REL=6.48
TRANSPOSONS_REL=10.2

REF_GENOME=ref/dmel-all-chromosome-r${FB_REL}.fasta
REF_TRANSPOSONS=ref/dmel-all-transposon-r${FB_REL}.fasta
TRANSPOSONS_DMEL=transposons/D_mel_transposon_sequence_set_v${TRANSPOSONS_REL}.fa


SX4_FOR=1_R1_001 2_R1_001 7_R1_001 8_R1_001
SX4_REV=1_R2_001 2_R2_001 7_R2_001 8_R2_001
SX4=${SX4_FOR} ${SX4_REV}



SX4_COV=${addprefix output/fast_coverage/, ${addsuffix .txt, ${SX4}}}
SX4_TE_MAP=${addprefix output/te_mapper/, ${addsuffix /te_mapper_output.json, 1_R1_001}}
MISC_1=output/te_mapper/1_R1_001/split_reads_for_tgt.csv

TARGETS=${SX4_COV} ${SX4_TE_MAP} ${MISC_1}

all: ${TARGETS}

clean:
	rm -rf output/*

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
output/te_mapper/1_R1_001/split_reads_for_tgt.csv: scripts/get_table1_split_reads.py output/te_mapper/1_R1_001/te_mapper_output.json
	python3 $<