#!/bin/bash

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --all_samples)
            all_samples=true
            shift
            ;;
        --include_samples)
            include_samples="$2"
            shift 2
            ;;
        --skip_samples)
            skip_samples="$2"
            shift 2
            ;;
        --rna)
            rna=true
            shift
            ;;
        --atac)
            atac=true
            shift
            ;;
        --souporcell)
            souporcell=true
            shift
            ;;
        --assign_geno)
            assign_geno=true
            shift
            ;;
        --overwrite_rna)
            overwrite_rna=true
            shift
            ;;
        -i|--input)
            input_directory="$2"
            shift 2
            ;;
        -o|--output)
            output_directory="$2"
            shift 2
            ;;
        -d|--dd_file)
            dd_file="$2"
            shift 2
            ;;
        -r|--ref_genome)
            ref_genome="$2"
            shift 2
            ;;
        -h|--help)
            help=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Display help message
if [[ $help == true ]]; then
    echo "Usage: sh multiome_qc_workflow.sh [OPTIONS]"
    echo "Example1: sh multiome_qc_workflow.sh --all_samples --rna --atac -i 'input_directory' -o 'output_directory'"
    echo "Example2: sh multiome_qc_workflow.sh --include_samples sample1,sample2,sample3 --rna -i 'input_directory' -o 'output_directory'"
    echo "Options:"
    echo "  --all_samples      Process all samples"
    echo "  --include_samples  Include specific samples (comma-separated)"
    echo "  --skip_samples     Skip specific samples (comma-separated)"
    echo "  --rna              Process RNA data"
    echo "  --atac             Process ATAC data"
    echo "  --souporcell       Run souporcell doublet detection on RNA data"
    echo "  --assign_geno      Correlate cluster to donor reference SNP genotypes"
    echo "  --overwrite_rna    Overwrite RNA preprocessing results."
    echo "  -i, --input        Input directory (required)"
    echo "  -o, --output       Output directory (required)"
    echo "  -d, --dd_file      Input csv file for souporcell double detection. Columns = sample_name,n_genotypes,vcf_file. (required to run souporcell)"
    echo "  -r, --ref_genome   Refernce genome directory for doublet detection"
    exit 0
fi

# Check for mandatory arguments
if [[ !($rna || $atac || $souporcell || $assign_geno) ]]; then
    echo "Error: Specify either RNA, ATAC, souporcell, or assign_geno sub-workflow. Refer to the --help function."
    exit 1
elif [[ $assign_geno != true && ( -z $input_directory || -z $output_directory ) ]]; then
    echo "Error: Input and output directories are required. Refer to the --help function."
    exit 1
fi

if [[ $assign_geno || $souporcell ]]; then
    if [[ -z $output_directory || -z $dd_file ]]; then
        echo "Error: Input CSV file and output folder are required. Refer to the --help function."
        exit 1
    fi
fi

# Function to submit LSF job
submit_job() {
    sample_name="$1"
    script_type="$2"
    input_directory="$3"
    output_directory="$4"
    
    # Find cellranger output directory
    if [[ $rna == true || $atac == true || $souporcell == true ]]; then
        if [[ -d "${input_directory}/${sample_name}/raw_feature_bc_matrix" ]]; then
            cr_outs="${input_directory}/${sample_name}"
        elif [[ -d "${input_directory}/${sample_name}/outs/raw_feature_bc_matrix" ]]; then
            cr_outs="${input_directory}/${sample_name}/outs"
        else
            echo "Error: raw_feature_bc_matrix directory not found for sample $sample_name"
            return 1
        fi
    fi
    
    # Create sample output directory
    if [ ! -d "${output_directory}/${sample_name}" ]; then
        mkdir -p "${output_directory}/${sample_name}"
        echo "Creating: ${output_directory}/${sample_name}"
    fi
    
    #-------------------------------------------------#
                 # PRE-PROCESS RNA DATA
    #-------------------------------------------------#
    if [[ $script_type == "RNA" ]]; then
        echo "Submitting job for $sample_name (type: $script_type)"
        
        overwrite_rna="$5"
        if [[ "$overwrite_rna" == true ]]; then
            overwrite_flag="--overwrite"
        else
            overwrite_flag=""
        fi
    
        # Create a temporary LSF script file for RNA
        temp_script=$(mktemp)
        cat > $temp_script <<EOF
#!/bin/bash
#BSUB -J "rna_qc_$sample_name"
#BSUB -q gpuqueue
#BSUB -W 5:00
#BSUB -n 2
#BSUB -R "rusage[mem=100]"
#BSUB -e "$output_directory/$sample_name/%J.err"
#BSUB -o "$output_directory/$sample_name/%J.out"

source ~/.bashrc
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/
conda run -n multiome_winner_test3 python rna_init.py -s "$sample_name" -i "$cr_outs" -o "$output_directory" $overwrite_flag
EOF

        # Submit the job using the temporary script file
        bsub < $temp_script

        # Remove the temporary script file
        rm $temp_script

    #-------------------------------------------------#
                # PRE-PROCESS ATAC DATA
    #-------------------------------------------------#
    elif [[ $script_type == "ATAC" ]]; then
        echo "Submitting job for $sample_name (type: $script_type)"
    
        # Create a temporary LSF script file for ATAC
        temp_script=$(mktemp)
        cat > $temp_script <<EOF
#!/bin/bash
#BSUB -J "atac_qc_$sample_name"
#BSUB -q gpuqueue
#BSUB -W 5:00
#BSUB -n 2
#BSUB -R "rusage[mem=100]"
#BSUB -e "$output_directory/$sample_name/%J.err"
#BSUB -o "$output_directory/$sample_name/%J.out"

source ~/.bashrc
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/
conda run -n multiome_winner_test3 Rscript atac_init.R "$sample_name" "$cr_outs" "$output_directory"
EOF

        # Submit the job using the temporary script file
        bsub < $temp_script

        # Remove the temporary script file
        rm $temp_script

    #-------------------------------------------------#
             # SOUPORCELL DOUBLET DETECTION
    #-------------------------------------------------#
    elif [[ $script_type == "souporcell" ]]; then
        echo "Submitting job for $sample_name (type: $script_type)"
        assign_geno=${5:-false}
        ref_genome="$6"
        n_genotypes="$7"
        vcf_file="$8"
        
        soup_out="$output_directory/$sample_name/RNAsoup_out"
        if [ ! -d "$soup_out" ]; then
            mkdir -p "$soup_out"
            echo "Creating: $soup_out"
        fi
        
        # Create a temporary LSF script file for ATAC
        temp_script=$(mktemp)
        cat > $temp_script <<EOF
#!/bin/bash
#BSUB -J "souporcell_$sample_name"
#BSUB -q gpuqueue
#BSUB -W 24:00
#BSUB -n 2
#BSUB -R "rusage[mem=100]"
#BSUB -e "$output_directory/$sample_name/%J.err"
#BSUB -o "$output_directory/$sample_name/%J.out"

source ~/.bashrc

module load singularity/3.7.1
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/souporcell
     
# ilc pre
singularity exec --bind "$cr_outs":/data,"$ref_genome":/ref,"$soup_out":/out souporcell_latest.sif souporcell_pipeline.py \
    -i /data/gex_possorted_bam.bam \
    -b /data/filtered_feature_bc_matrix/barcodes.tsv.gz \
    -f /ref/fasta/genome.fa \
    -t 8 \
    -o /out/ \
    -k "$n_genotypes"
    
singularity exec --bind "$soup_out":/out Demuxafy.sif bash souporcell_summary.sh /out/clusters.tsv > "$soup_out"/souporcell_summary.tsv

if [[ $assign_geno == true ]]; then
    singularity exec --bind "$soup_out":/out,"$vcf_file":/mnt/genotypes.vcf Demuxafy.sif Assign_Indiv_by_Geno.R \
        -r /mnt/genotypes.vcf \
        -c /out/cluster_genotypes.vcf \
        -o /out
fi
EOF

        # Submit the job using the temporary script file
        bsub < $temp_script

        # Remove the temporary script file
        rm $temp_script
        
    
    #-------------------------------------------------#
             # DEMUXAFY GENOTYPE ASSIGNMENT
    #-------------------------------------------------#
    elif [[ $script_type == "assign_geno" ]]; then
        echo "Submitting job for $sample_name (type: $script_type)"
        vcf_file="$5"
        
        soup_out="$output_directory/$sample_name/RNAsoup_out"
        if [[ -d "$soup_out/cluster_genotypes.vcf" ]]; then
            echo "Error: $soup_out/cluster_genotypes.vcf not found"
            return 1
        elif [[ -d "$vcf_file" ]]; then
            echo "Error: $vcf_file not found"
            return 1
        fi
        
        # Create a temporary LSF script file for ATAC
        temp_script=$(mktemp)
        cat > $temp_script <<EOF
#!/bin/bash
#BSUB -J "assign_geno_$sample_name"
#BSUB -q gpuqueue
#BSUB -W 24:00
#BSUB -n 2
#BSUB -R "rusage[mem=100]"
#BSUB -e "$output_directory/$sample_name/%J.err"
#BSUB -o "$output_directory/$sample_name/%J.out"

source ~/.bashrc

module load singularity/3.7.1
cd /data/niecr/ramaiar/multiome_lucyrepo/analysis/souporcell
     
singularity exec --bind "$soup_out":/out,"$vcf_file":/mnt/genotypes.vcf Demuxafy.sif Assign_Indiv_by_Geno.R \
    -r /mnt/genotypes.vcf \
    -c /out/cluster_genotypes.vcf \
    -o /out
EOF

        # Submit the job using the temporary script file
        bsub < $temp_script

        # Remove the temporary script file
        rm $temp_script

    fi
}

# Read all sample names in the input directory
sample_names=($(ls "$input_directory"))

# Process samples based on the command-line arguments
if [[ $all_samples == true ]]; then
    for sample_name in "${sample_names[@]}"; do
        if [[ $rna == true ]]; then
            submit_job "$sample_name" "RNA" "$input_directory" "$output_directory" "$overwrite_rna"
        fi
        if [[ $atac == true ]]; then
            submit_job "$sample_name" "ATAC" "$input_directory" "$output_directory"
        fi
    done
elif [[ -n $include_samples ]]; then
    include_samples_array=($(echo "$include_samples" | tr ',' '\n'))
    for sample_name in "${include_samples_array[@]}"; do
        if [[ " ${sample_names[@]} " =~ " $sample_name " ]]; then
            if [[ $rna == true ]]; then
                submit_job "$sample_name" "RNA" "$input_directory" "$output_directory" "$overwrite_rna"
            fi
            if [[ $atac == true ]]; then
                submit_job "$sample_name" "ATAC" "$input_directory" "$output_directory"
            fi
        else
            echo "Sample $sample_name not found in the input directory."
        fi
    done
elif [[ -n $skip_samples ]]; then
    skip_samples_array=($(echo "$skip_samples" | tr ',' '\n'))
    for sample_name in "${sample_names[@]}"; do
        if [[ ! " ${skip_samples_array[@]} " =~ " $sample_name " ]]; then
            if [[ $rna == true ]]; then
                submit_job "$sample_name" "RNA" "$input_directory" "$output_directory" "$overwrite_rna"
            fi
            if [[ $atac == true ]]; then
                submit_job "$sample_name" "ATAC" "$input_directory" "$output_directory"
            fi
        fi
    done
fi

if [[ $souporcell == true || $assign_geno == true ]]; then
    # Path to your CSV file
    input_csv_file=$dd_file

    # Read the CSV file line by line
    while IFS=, read -r sample_name n_genotypes vcf_file; do
        # Skip the header line
        if [[ $sample_name == "sample_name" ]]; then
            continue
        fi

        echo "Processing sample: $sample_name"
        echo "Number of genotypes: $n_genotypes"
        echo "VCF files: $vcf_file"

        if [[ $souporcell == true ]]; then
            submit_job "$sample_name" "souporcell" "$input_directory" "$output_directory" "$assign_geno" "$ref_genome" "$n_genotypes" "$vcf_file"
        elif [[ $assign_geno == true ]]; then
            submit_job "$sample_name" "assign_geno" "$input_directory" "$output_directory" "$vcf_file"
        fi
    done < "$input_csv_file"
fi
