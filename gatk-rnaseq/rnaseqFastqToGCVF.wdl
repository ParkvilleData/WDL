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
  String genomeDir
  String opossumLocation
  File dbSnpVcf
  File dbSnpVcfIndex
  File knownVcfs
  File knownVcfsIndices

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
	Int convertSamToBamRunThreads
        Int convertSamToBamRunMinutes
        Int convertSamToBamRunMem

                input:
                   alignmentSam=starAlignment.outputSam,
                   sampleName=sample[0]
        }

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
           picardInputBam=convertSamToBam.outputBam,
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


	call sortBam {
                input:
                        bam2sort=splitNCigarReads.splitCigarsBamOutput,
                        sampleName=sample[0]
        }


	call baseReCalibrator1 {
                input:
                        gatkLocation=gatkLocation,
                        sortedBam=sortBam.sortedBam,
                        dbsnp=dbSnpVcf,
                        goldStandard=knownVcfs,
                        sampleName=sample[0],
                        ref_fasta=refFasta,
                        ref_fasta_index=refIndex,
                        ref_dict=refDict,
                        outputSortedBamIndex=sortBam.outputSortedBam
        }

        call baseReCalibrator2 {
                input:
                        gatkLocation=gatkLocation,
                        sortedBam=sortBam.sortedBam,
                        dbsnp=dbSnpVcf,
                        goldStandard=knownVcfs,
                        sampleName=sample[0],
                        ref_fasta=refFasta,
                        ref_fasta_index=refIndex,
                        ref_dict=refDict,
                        outputSortedBamIndex=sortBam.outputSortedBam,
                        calibratedGrp=baseReCalibrator1.calibratedFile1
        }

        call generatePlots {
                input:
                        gatkLocation=gatkLocation,
                        ref_fasta=refFasta,
                        ref_fasta_index=refIndex,
                        sampleName=sample[0],
                        ref_dict=refDict,
                        calibratedFile1=baseReCalibrator1.calibratedFile1,
                        calibratedFile2=baseReCalibrator2.calibratedFile2
        }

        call printReads {
                input:
                        gatkLocation=gatkLocation,
                        ref_fasta=refFasta,
                        ref_fasta_index=refIndex,
                        ref_dict=refDict,
                        sortedBam=splitNCigarReads.splitCigarsBamOutput,
                        calibratedFile1=baseReCalibrator1.calibratedFile1,
                        sortedBamIndex=createBamIndex.splitCigarsBamIndexOutput,
                        sampleName=sample[0]
        }

        call indexCalibratedBam {
                input:
                        sampleName=sample[0],
			refFasta=refFasta,
                        calBam=printReads.recalibratedReadsBam
        }

	call opossum {
		input:
		 calibratedBam=indexCalibratedBam.calibratedBam,
		 opossumLocation=opossumLocation,
		 sampleName=sample[0]

	}

	call indexOpossumBam {
		input:
		 sampleName=sample[0],
		 opossumBam=opossum.opossumBam
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
        bamFile=indexOpossumBam.calibratedBam, 
        bamIndex=indexOpossumBam.calibratedBamIndex
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
#end workflow calls

#start tasks
 task starAlignment {
        Int starAlignmentRunThreads
        Int starAlignmentRunMinutes
        Int starAlignmentRunMem
        File inputFastqRead1
        File inputFastqRead2
	String genomeDir
	String sampleName

        command {
                module load STAR

                STAR --genomeDir ${genomeDir} \
                        --runThreadN ${starAlignmentRunThreads} \
                        --readFilesIn ${inputFastqRead1} ${inputFastqRead2} \
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
        Int convertSamToBamRunThreads
        Int convertSamToBamRunMinutes
        Int convertSamToBamRunMem
        String sampleName

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

       task sortBam {
                File bam2sort
                String sampleName
		Int sortBamMem
		Int sortBamThreads
		Int sortBamRunMinutes

                command {
			module load SAMtools
                        samtools sort ${bam2sort} -o ${sampleName}.sorted.bam
                        samtools index ${sampleName}.sorted.bam
                }
		runtime {
			runtime_minutes: '${sortBamRunMinutes}'
			cpus: '${sortBamThreads}'
			mem: '${sortBamMem}'
		}
                output {
                        File sortedBam="${sampleName}.sorted.bam"
                        File outputSortedBam="${sampleName}.sorted.bam.bai"
                }
        }



 #first round of GATK BaseCalibrator
        task baseReCalibrator1 {
                File gatkLocation
                File sortedBam
                File dbsnp
                File goldStandard
                File ref_fasta
                File ref_fasta_index
                File ref_dict
                File outputSortedBamIndex
                String sampleName
		Int baseRecal1RunMinutes
		Int baseRecal1Threads
		Int baseRecal1Mem

                command {
			module load Java
                       java -jar ${gatkLocation} -T BaseRecalibrator \
                                -R ${ref_fasta} \
                                -I ${sortedBam} \
                                -o ${sampleName}.calibrated.grp \
                                -knownSites ${dbsnp} \
                                -knownSites ${goldStandard}
                }
		runtime {
			runtime_minutes: '${baseRecal1RunMinutes}'
			cpus: '${baseRecal1Threads}'
			mem: '${baseRecal1Mem}'
		}
                output {
                        File calibratedFile1="${sampleName}.calibrated.grp"
                }
        }

        #second pass of recalibrator
        task baseReCalibrator2 {
                File gatkLocation
                File sortedBam
                File dbsnp
                File goldStandard
                File ref_fasta
                File ref_fasta_index
                File ref_dict
                File outputSortedBamIndex
                File calibratedGrp
                String sampleName
		Int baseRecal2RunMinutes
		Int baseRecal2Threads
		Int baseRecal2Mem

                command {
			module load Java

                        java -jar ${gatkLocation} -T BaseRecalibrator \
                                -R ${ref_fasta} \
                                -I ${sortedBam} \
                                -knownSites ${dbsnp} \
                                -knownSites ${goldStandard} \
                                -BQSR ${calibratedGrp} \
                                -o ${sampleName}.post.calibrated.grp
                }
		runtime {
			runtime_minutes: '${baseRecal2RunMinutes}'
			cpus: '${baseRecal2Threads}'
			mem: '${baseRecal2Mem}'
		}
                output {
                        File calibratedFile2="${sampleName}.post.calibrated.grp"
                }
        }
 #generate before and after plots
        task generatePlots {
                File gatkLocation
                File ref_fasta
                File ref_fasta_index
                File ref_dict
                File calibratedFile1
                File calibratedFile2
                String sampleName
		Int generatePlotsRunMinutes
		Int generatePlotsThreads
		Int generatePlotsMem

                command {
			module load R/3.2.3-GCC-4.9.3
			module load Java
                        java -jar ${gatkLocation} -T AnalyzeCovariates \
                                -R ${ref_fasta} \
                                -before ${calibratedFile1} \
                                -after ${calibratedFile2} \
                                -plots ${sampleName}.plots.pdf
                }
		runtime {
			runtime_minutes: '${generatePlotsRunMinutes}'
			cpus: '${generatePlotsThreads}'
			mem: '${generatePlotsMem}'
		}
                output {
                        File calibratedPlots="${sampleName}.plots.pdf"
                }
        }


        #print reads
        task printReads {
                File gatkLocation
                File ref_fasta
                File ref_fasta_index
                File ref_dict
                File sortedBam
                File calibratedFile1
                File sortedBamIndex
                String sampleName
		Int printReadsRunMinutes
		Int printReadsThreads
		Int printReadsMem

                command {
			module load Java

                        java -jar ${gatkLocation} -T PrintReads \
                                -R ${ref_fasta} \
                                -I ${sortedBam} \
                                -BQSR ${calibratedFile1} \
                                -o ${sampleName}.recal_reads.bam
                }
		runtime {
			runtime_minutes: '${printReadsRunMinutes}'
			cpus: '${printReadsThreads}'
			mem: '${printReadsMem}'
		}
                output {
                       File recalibratedReadsBam="${sampleName}.recal_reads.bam"
                }

        }

        #index recalibrated bam
        task indexCalibratedBam {
                File calBam
		File refFasta
                String sampleName
		Int indexCalibratedBamRunMinutes
		Int indexCalibratedBamThreads
		Int indexCalibratedBamMem

                command {
			module load SAMtools

			samtools calmd ${calBam} ${refFasta} -b > ${sampleName}.temp.bam
                        samtools sort ${sampleName}.temp.bam -o ${sampleName}.recal_reads.sorted.bam
                        samtools index ${sampleName}.recal_reads.sorted.bam
                }
		runtime {
			runtime_minutes: '${indexCalibratedBamRunMinutes}'
			cpus: '${indexCalibratedBamThreads}'
			mem: '${indexCalibratedBamMem}'
		}
                output {
                        File calibratedBam="${sampleName}.recal_reads.sorted.bam"
                        File calibratedBamIndex="${sampleName}.recal_reads.sorted.bam.bai"
                }

	}

  task opossum {
	File calibratedBam
	Int opossumRunMinutes
	Int opossumThreads
	Int opossumMem
	String opossumLocation
	String sampleName

	command {
                        module load Python/2.7.12-GCC-4.9.3

			python ${opossumLocation} --BamFile=${calibratedBam} --OutFile=${sampleName}.opossum.bam
                }
                runtime {
                        runtime_minutes: '${opossumRunMinutes}'
                        cpus: '${opossumThreads}'
                        mem: '${opossumMem}'
                }
                output {
                        File opossumBam="${sampleName}.opossum.bam"
                }	
  }

  task indexOpossumBam {
	File opossumBam
	Int indexOpossumBamRunMinutes
	Int indexOpossumBamThreads
	Int indexOpossumBamMem
	String sampleName

	
	command {
                        module load SAMtools

                        samtools sort ${opossumBam} -o ${sampleName}.opossum.sorted.bam
                        samtools index ${sampleName}.opossum.sorted.bam
                }
                runtime {
                        runtime_minutes: '${indexOpossumBamRunMinutes}'
                        cpus: '${indexOpossumBamThreads}'
                        mem: '${indexOpossumBamMem}'
                }
                output {
                        File calibratedBam="${sampleName}.opossum.sorted.bam"
                        File calibratedBamIndex="${sampleName}.opossum.sorted.bam.bai"
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
	module load Java

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
	Int variantFilterRunMinutes
	Int variantFilterThreads
	Int variantFilterMem


	command {
	module load Java

		java -jar  ${gatk_path} \
		        -T VariantFiltration \
			-R ${ref_fasta} \
			-V ${input_vcf} \
			-window 35 \
			-cluster 3 \
			--filterName "FS" \
			-filter "FS > 30.0" \
			--filterName "QD" \
			-filter "QD < 2.0" \
			-o ${base_name}.filtered.vcf
	}

	output {
    	File output_vcf = "${base_name}.filtered.vcf"
    	File output_vcf_index = "${base_name}.filtered.vcf.idx"
	}

	runtime {
	   runtime_minutes: '${variantFilterRunMinutes}'
           cpus: '${variantFilterThreads}'
           mem: '${variantFilterMem}'
	}
}
