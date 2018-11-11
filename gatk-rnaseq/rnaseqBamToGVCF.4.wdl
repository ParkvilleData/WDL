## rnaseq GATK best practices pipeline
## Bobbie Shaban
## 09/04/2018
## Melbourne Integrative genomics
##
## This GATK pipeline will take a number of bam files
## and convert them to GCVF files and the output
## it will employ scatter gather and hpc
## instructions to make are here: https://gatkforums.broadinstitute.org/wdl/discussion/6716/scatter-gather-parallelism
## and here: https://gatkforums.broadinstitute.org/wdl/discussion/7614/4-howto-use-scatter-gather-to-joint-call-genotypes

workflow rnaseqBamToGCVF {

  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  File refFasta
  File refDict
  File refIndex
  String picardLocation
  String gatkLocation
  File dbSnpVcf
  File dbSnpVcfIndex
  File knownVcfs
  File knownVcfsIndices

  scatter (sample in inputSamples) {
    call picard {
       String typeARRG
       String sortOrder
       String readGroupID
       String readGroupLibrary
       String readGroupPlatform
       String readGroupPlatformBarcode
       Int picardRunMinutes
       Int picardThreads

       input:
           picardLocation=picardLocation,
           picardInputBam=sample[1],
           sampleName=sample[0]
    }

    call picardMarkDuplicates {
           String typeMD
           String createIndex
           String validationStringency
           String outputMetrics
 	   String picardMarkDuplicatesRunMinutes
	   String picardMarkDuplicatesThreads

            input:
                picardLocation=picardLocation,
                picardMDInputBam=picard.picardOutputBam,
                sampleName=sample[0]
    }

    call createPicardBamIndex {
   	  Int createPicardBamIndexRunMinutes
	  Int createPicardBamIndexThreads
	  Int createPicardBamIndexMem

            input:
                picardMDBamToBeIndexed=picardMarkDuplicates.picardDeduppedBam
    }

    call createRefIndex{
	  Int createRefIndexRunMinutes
	  Int createRefIndexRunThreads
	  Int createRefIndexMem
           input:
                refFasta=refFasta
    }

    call splitNCigarReads{
            String splitCigars
            Int RF
            Int RMQF
            String RMQT
            String U
	    Int splitNCigarReadsRunMinutes	
	    Int splitNCigarReadsThreads
	    Int splitNCigarReadsMem

            input:
                refFasta=refFasta,
                refFastaIndex=createRefIndex.refFastaIndex,
                refDictionary=refDict,
                sampleName=sample[0],
                gatkLocation=gatkLocation,
                splitCigarsInputBam=picardMarkDuplicates.picardDeduppedBam,
                splitCigarsInputBamIndex=createPicardBamIndex.mDBamIndex
        }

    call createBamIndex{
	    Int createBamIndexRunMinutes
	    Int createBamIndexThreads
	    Int createBamIndexMem

            input:
                bamToBeIndexed=splitNCigarReads.splitCigarsBamOutput
    }

    call BaseRecalibrator {
	   Int BaseRecalRunMinutes
	   Int BaseRecalRunThreads
	   Int BaseRecalMem

          input:
                input_bam = splitNCigarReads.splitCigarsBamOutput,
                input_bam_index = createBamIndex.splitCigarsBamIndexOutput,
                recal_output_file = sample[0] + ".recal_data.csv",
                dbSNP_vcf = dbSnpVcf,
                dbSNP_vcf_index = dbSnpVcfIndex,
                known_indels_sites_VCFs = knownVcfs,
                known_indels_sites_indices = knownVcfsIndices,
                ref_dict = refDict,
                ref_fasta = refFasta,
                ref_fasta_index = refIndex,
                gatk_path = gatkLocation
    }

    call ApplyBQSR {
	Int ApplyBQSRRunMinutes
	Int ApplyBQSRThreads
	Int ApplyBQSRMem

		input:
			input_bam = splitNCigarReads.splitCigarsBamOutput,
			input_bam_index = createBamIndex.splitCigarsBamIndexOutput,
			base_name = sample[0] + ".aligned.duplicates_marked.recalibrated",
			ref_fasta = refFasta,
			ref_fasta_index = refIndex,
			ref_dict = refDict,
			recalibration_report = BaseRecalibrator.recalibration_report,
			gatk_path = gatkLocation
   }

   call bqsrBamIndex {
	Int bqsrBamIndexRunMinutes
	Int bqsrBamIndexThreads
	Int bqsrBamIndexMem
 
	input: 
		bamToBeIndexed = ApplyBQSR.output_bam
   }

    call HaplotypeCallerERC {
	Int haplotypeCallerRunMinutes
	Int haplotypeCallerThreads
	Int haplotypeCallerMem

      input: GATK=gatkLocation, 
        RefFasta=refFasta, 
	RefIndex=refIndex,
        RefDict=refDict, 
        sampleName=sample[0],
        bamFile=ApplyBQSR.output_bam, 
        bamIndex=bqsrBamIndex.bqsrBamOutput
    }
  }

  call combineGVCFs {
   Int combineRunMinutes
   Int combineRunThreads
   Int combineRunMem

   input: GATK=gatkLocation,
	RefFasta=refFasta,
	RefIndex=refIndex,
	RefDict=refDict,
	GVCFs=HaplotypeCallerERC.GVCF,
	sampleName="combinedGVCFs"	
  }

  call GenotypeGVCFs {
    Int genotypeRunMinutes
    Int genotypeThreads
    Int genotypeMem

    input: GATK=gatkLocation, 
      RefFasta=refFasta, 
      RefIndex=refIndex, 
      RefDict=refDict, 
      sampleName="CEUtrio", 
      combinedVCF=combineGVCFs.combinedOutput
  }

  call VariantFiltration {
	Int variantFilterRunMinutes
	Int variantFilterThreads
	Int variantFilterMem

	input:
		input_vcf = GenotypeGVCFs.rawVCF,
		input_vcf_index = GenotypeGVCFs.rawVCFidx,
		base_name = "CEUtrio",
		ref_fasta = refFasta,
		ref_fasta_index = refIndex,
		ref_dict = refDict,
		gatk_path=gatkLocation
	}

}

 task picard {
        String typeARRG
        String sampleName
        String sortOrder
        String readGroupID
        String readGroupLibrary
        String readGroupPlatform
        String readGroupPlatformBarcode
        String picardLocation
	Int picardRunMinutes
	Int picardThreads
	Int picardMem
        File picardInputBam

        command {
		module load Java

                java -jar ${picardLocation} ${typeARRG} \
                        I=${picardInputBam} \
                        O=${sampleName}Aligned.bam \
                        SO=${sortOrder} \
                        RGID=${readGroupID} \
                        RGLB=${readGroupLibrary} \
                        RGPL=${readGroupPlatform} \
                        RGPU=${readGroupPlatformBarcode} \
                        RGSM=${sampleName}
        }
	runtime {
		runtime_minutes: '${picardRunMinutes}'
		cpus: '${picardThreads}'
		mem: '${picardMem}'
	}
        output {
                File picardOutputBam="${sampleName}Aligned.bam"
        }
 }

 task picardMarkDuplicates {
        String typeMD
        String sampleName
        String picardLocation
        String createIndex
        String validationStringency
        String outputMetrics
        File picardMDInputBam
	Int picardMarkDuplicatesRunMinutes
	Int picardMarkDuplicatesThreads
	Int picardMarkDuplicatesMem

        command {
		module load Java
                java -jar ${picardLocation} ${typeMD} \
                     I=${picardMDInputBam} \
                     O="${sampleName}.dedupped.bam" \
                     CREATE_INDEX=${createIndex} \
                     VALIDATION_STRINGENCY=${validationStringency} \
                     M=${outputMetrics}
        }
	runtime {
                runtime_minutes: '${picardMarkDuplicatesRunMinutes}'
		cpus: '${picardMarkDuplicatesThreads}'
		mem: '${picardMarkDuplicatesMem}'
       }
        output {
                File picardDeduppedBam="${sampleName}.dedupped.bam"
        }
 }

 task createPicardBamIndex {
        File picardMDBamToBeIndexed
	Int createPicardBamIndexRunMinutes
	Int createPicardBamIndexThreads
	Int createPicardBamIndexMem

        command {
		module load SAMtools

                samtools index ${picardMDBamToBeIndexed}
        }
	runtime {
                runtime_minutes: '${createPicardBamIndexRunMinutes}'
		cpus: '${createPicardBamIndexThreads}'
		mem: '${createPicardBamIndexMem}'
        }
        output {
                File mDBamIndex=sub("${picardMDBamToBeIndexed}", "bam", "bai")
        }
 }

 task createRefIndex {
        File refFasta
	Int createRefIndexRunMinutes
	Int createRefIndexRunThreads
	Int createRefIndexMem

        command {
		module load SAMtools

                samtools faidx ${refFasta}
        }
	runtime {
                runtime_minutes: '${createRefIndexRunMinutes}'
		cpus: '${createRefIndexRunThreads}'
		mem: '${createRefIndexMem}'
        }
        output{
                File refFastaIndex = "${refFasta}.fai"
        }
 }

 task splitNCigarReads{
        String splitCigars
        String sampleName
        String gatkLocation
        String RF
        Int RMQF
        Int RMQT
        String U
        File refFasta
        File refFastaIndex
        File refDictionary
        File splitCigarsInputBam
        File splitCigarsInputBamIndex
	Int splitNCigarReadsRunMinutes
	Int splitNCigarReadsThreads
	Int splitNCigarReadsMem

        command {
		module load Java

                java -jar ${gatkLocation} \
                        -T ${splitCigars} \
                        -R ${refFasta} \
                        -I ${splitCigarsInputBam} \
                        -o ${sampleName}.splitCigars.bam \
                        -rf ${RF} \
                        -RMQF ${RMQF} \
                        -RMQT ${RMQT} \
                        --allow_potentially_misencoded_quality_scores \
                        -U ${U}
        }
	runtime {
                runtime_minutes: '${splitNCigarReadsRunMinutes}'
		cpus: '${splitNCigarReadsThreads}'
		mem: '${splitNCigarReadsMem}'
        }
        output {
                File splitCigarsBamOutput="${sampleName}.splitCigars.bam"
        }
 }

 task createBamIndex {
        File bamToBeIndexed
        Int createBamIndexRunMinutes 
	Int createBamIndexThreads
	Int createBamIndexMem

        command {
		module load SAMtools

                samtools index ${bamToBeIndexed}
        }
	runtime {
                runtime_minutes: '${createBamIndexRunMinutes}'
		cpus: '${createBamIndexThreads}'
		mem: '${createBamIndexMem}'
        }
        output {
                File splitCigarsBamIndexOutput=sub("${bamToBeIndexed}", "bam", "bai")
        }
 }

 task BaseRecalibrator {
    File input_bam
    File input_bam_index
    String recal_output_file
    File dbSNP_vcf
    File dbSNP_vcf_index
    File known_indels_sites_VCFs
    File known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String gatk_path
    Int BaseRecalRunMinutes
    Int BaseRecalThreads
    Int BaseRecalMem

    command {
	module load Java

        java -jar ${gatk_path} \
            -T BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -o ${recal_output_file} \
            -knownSites ${dbSNP_vcf} \
            -knownSites ${sep=" --known-sites " known_indels_sites_VCFs}
    }

    output {
        File recalibration_report = recal_output_file
    }

    runtime {
	runtime_minutes: '${BaseRecalRunMinutes}' 	
	cpus: '${BaseRecalThreads}'
	mem: '${BaseRecalMem}'
    }
}

task ApplyBQSR {
    File input_bam
    File input_bam_index
    String base_name
    File recalibration_report
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String gatk_path
    Int ApplyBQSRRunMinutes
    Int ApplyBQSRThreads
    Int ApplyBQSRMem

    command {
	module load Java

        java -jar ${gatk_path} \
            -T ApplyBQSR \
            --add-output-sam-program-record \
            -R ${ref_fasta} \
            -I ${input_bam} \
            -O ${base_name}.bam \
            --bqsr-recal-file ${recalibration_report}
    }

    output {
        File output_bam = "${base_name}.bam"
    }

    runtime {
	runtime_minutes: '${ApplyBQSRRunMinutes}'
        cpus: '${ApplyBQSRThreads}'
        mem: '${ApplyBQSRMem}'
    }
}

  task bqsrBamIndex {
        File bamToBeIndexed
        Int bqsrBamIndexRunMinutes
        Int bqsrBamIndexThreads
        Int bqsrBamIndexMem

        command {
                module load SAMtools
                samtools index ${bamToBeIndexed}
        }
        runtime {
                runtime_minutes: '${bqsrBamIndexRunMinutes}'
                cpus: '${bqsrBamIndexThreads}'
                mem: '${bqsrBamIndexMem}'
        }
        output {
                File bqsrBamOutput=sub("${bamToBeIndexed}", "bam", "bai")
        }
 }
 
  task HaplotypeCallerERC {

  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  Int haplotypeCallerThreads
  Int haplotypeCallerRunMinutes
  Int haplotypeCallerMem
  File bamFile
  File bamIndex

  command {
    module load Java

    java -jar ${GATK} \
        -T HaplotypeCaller \
	-ERC GVCF \
        -R ${RefFasta} \
        -I ${bamFile} \
        -o ${sampleName}_rawLikelihoods.g.vcf \
        -nct ${haplotypeCallerThreads}
  }
  runtime {
          runtime_minutes: '${haplotypeCallerRunMinutes}'
	  cpus: '${haplotypeCallerThreads}'
	  mem: '${haplotypeCallerMem}'	  
  }
  output {
   	  File GVCF = "${sampleName}_rawLikelihoods.g.vcf"
  }
 }

 task combineGVCFs {
  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  String sampleName
  Int combineRunMinutes
  Int combineRunThreads
  Int combineRunMem

  Array[File] GVCFs
 
  command {
	java -jar ${GATK} \
	   -T CombineGVCFs \
	   -R ${RefFasta} \
 	   --variant ${sep=" --variant " GVCFs} \
	   -o ${sampleName}.cohort.g.vcf		
  }
  runtime {
	runtime_minutes: '${combineRunMinutes}'
	cpus: '${combineRunThreads}'
	mem: '${combineRunMem}'
  }
  output {
     File combinedOutput = "${sampleName}.cohort.g.vcf"	
  }
 }


 task GenotypeGVCFs {

  File GATK
  File RefFasta
  File RefIndex
  File RefDict
  File combinedVCF
  String sampleName
  Int genotypeRunMinutes
  Int genotypeThreads
  Int genotypeMem

  command { 
    module load Java 

    java -Xmx100g -jar ${GATK} \
        -T GenotypeGVCFs \
        -R ${RefFasta} \
	-nt ${genotypeThreads} \
        -V ${combinedVCF} \
        -o ${sampleName}_rawVariants.vcf
  }
   runtime {
           runtime_minutes: '${genotypeRunMinutes}'
	   cpus: '${genotypeThreads}'
	   mem: '${genotypeMem}'
        }
  output {
    File rawVCF = "${sampleName}_rawVariants.vcf"
    File rawVCFidx = "${sampleName}_rawVariants.vcf.idx"
  }
}

task VariantFiltration {
	File input_vcf
	File input_vcf_index
	String base_name
  	File ref_dict
  	File ref_fasta
  	File ref_fasta_index
    	String gatk_path


	command {
		 ${gatk_path} \
		    VariantFiltration \
			--R ${ref_fasta} \
			--V ${input_vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" \
			--filter "FS > 30.0" \
			--filter-name "QD" \
			--filter "QD < 2.0" \
			-O ${base_name}.filtered.vcf
	}

	output {
    	File output_vcf = "${base_name}"
    	File output_vcf_index = "${base_name}.tbi"
	}

	runtime {
	   runtime_minutes: '${variantFilterRunMinutes}'
           cpus: '${variantFilterThreads}'
           mem: '${variantFilterMem}'
	}
}
