.PHONY: all clean



FB_REL=6.48

REF_GENOME=ref/dmel-all-chromosome-r${FB_REL}.fasta
REF_TRANSPOSONS=ref/dmel-all-transposon-r${FB_REL}.fasta



SX4_FOR=1_R1_001 2_R1_001 7_R1_001 8_R1_001
SX4_REV=1_R2_001 2_R2_001 7_R2_001 8_R2_001
SX4=${SX4_FOR} ${SX4_REV}



SX4_COV=${addprefix output/fast_coverage/, ${addsuffix .txt, ${SX4}}}

TARGETS=${SX4_COV}

all: ${TARGETS}

clean:
	rm -rf output/*

output/fast_coverage/%.txt:
	@mkdir -p ${dir $@}
	cargo run --bin calculate_stats --release -- fast_coverage -q reads/${basename ${notdir $@}}.fastq -r ${REF_GENOME} -o $@