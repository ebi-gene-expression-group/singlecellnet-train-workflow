#!/usr/bin/env nextflow  

// set up channels
TRAIN_DIR = Channel.fromPath(params.training_10x_dir)
TRAIN_METADATA = Channel.fromPath(params.training_metadata)

// if necessary, down-sample cells to avoid memory issues 
process downsample_cells {
    conda "${baseDir}/envs/label_analysis.yaml"

    memory { 32.GB * task.attempt }
    maxRetries 5
    errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }

    input:
        file(expression_data) from TRAIN_DIR
        file(training_metadata) from TRAIN_METADATA
        
    output:
        file("expr_data_downsampled") into TRAIN_DIR_DOWNSAMPLED
        file("metadata_filtered.tsv") into TRAIN_METADATA_DOWNSAMPLED

    """
    set +e
    downsample_cells.R\
        --expression-data ${expression_data}\
        --metadata ${training_metadata}\
        --exclusions ${params.exclusions}\
        --cell-id-field "${params.cell_id_col}"\
        --cell-type-field "${params.cell_types_col}"\
        --output-dir expr_data_downsampled\
        --metadata-upd metadata_filtered.tsv
    if [ \$? -eq  2 ];
    then
        cp -P ${expression_data} expr_data_downsampled
        cp -P ${training_metadata} metadata_filtered.tsv
        exit 0
    fi
    """
}

// produce sce object for training dataset
process create_training_sce {
    conda "${baseDir}/envs/dropletutils.yaml"

    memory { 32.GB * task.attempt }
    maxRetries 5
    errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }
    
    input:
        file(train_metadata) from TRAIN_METADATA_DOWNSAMPLED
        file(train_dir) from TRAIN_DIR_DOWNSAMPLED

    output:
        file("training_sce.rds") into TRAINING_SCE

    """
    dropletutils-read-10x-counts.R\
                --samples ${train_dir}\
                --col-names ${params.col_names}\
                --metadata-files ${train_metadata}\
                --cell-id-column "${params.cell_id_col}"\
                --metadata-columns "${params.cell_id_col}","${params.cell_types_col}"\
                --output-object-file training_sce.rds
    """ 
}

process train_scn_classifer {

    publishDir "${params.results_dir}"
    conda "${baseDir}/envs/singlecellnet.yaml"

    memory { 32.GB * task.attempt }
    maxRetries 5
    errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }

    input:
        file(sce_object) from TRAINING_SCE

    output:
        file("trained_classifer.rds") into TRAINED_CLASSIFIER

    """
    scn-train-model.R\
            --input-object ${sce_object}\
            --cell-type-col ${params.cell_types_col}\
            --cell-barcode-col ${params.cell_id_col}\
            --n-top-genes ${params.n_top_genes}\
            --n-top-gene-pairs ${params.n_top_gene_pairs}\
            --n-rand ${params.n_rand}\
            --n-trees ${params.n_trees}\
            --stratify ${params.stratify}\
            --weighted-down-threshold ${params.weighted_down_threshold}\
            --transprop-factor ${params.transprop_factor}\
            --weighted-down-total ${params.wieghted_down_total}\
            --dataset-id ${params.training_dataset_id}\
            --output-path trained_classifer.rds
    """
}

