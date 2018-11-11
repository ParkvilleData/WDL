workflow rnaseqFastqToBam {

  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File refFasta
  File refDict
  File refIndex
  String picardLocation
  String gatkLocation
  String genomeDir

  scatter (sample in inputSamples) {

	call starAlignment {
	 Int starAlignmentRunThreads
	 Int starAlignmentRunMinutes
	 Int starAlignmentRunMem

		input:
		   genomeDir=genomeDir,
		   inputFastqRead1=sample[1],
		   inputFastqRead2=sample[2],
		   sampleName=sample[0]
	}

	call convertSamToBam {
		input:
		   alignmentSam=starAlignment.outputSam,
		   sampleName=sample[0]
	}
  }

} #end of workflow

  task starAlignment {
	Int starAlignmentRunThreads
	Int starAlignmentRunMinutes
	Int starAlignmentRunMem
	String genomeDir
	String sampleName
	File inputFastqRead1
	File inputFastqRead2

	command {
                module load STAR

		STAR --genomeDir ${genomeDir} \
			--runThreadN ${starAlignmentRunThreads} \
			--readFilesIn ${inputFastqRead1} $inputFastqRead2 \
			--outFileNamePrefix ${sampleName}
        }
        runtime {
		runtime_minutes: '${starAlignmentRunMinutes}'
		cpus: '${starAlignmentRunThreads}'
		mem: '${starAlignmentRunMem}'
	}
	output {
		File outputSam="${sampleName}Aligned.out.sam"
	}
  }

  task convertSamToBam {
	File alignmentSam
	String sampleName
	Int convertSamToBamRunThreads
	Int convertSamToBamRunMinutes
	Int convertSamToBamRunMem

	command {
		module load SAMtools
	
		samtools view -bS ${alignmentSam} > ${sampleName}.aligned.bam
	}
	runtime {
		runtime_minutes: '${convertSamToBamRunMinutes}'
		cpus: '${convertSamToBamRunThreads}'
		mem: '${convertSamToBamRunMem}'
	}
	output {
		File outputBam = "${sampleName}.aligned.bam"
	}
  }
