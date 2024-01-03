#!/usr/bin/env python3


import argparse
import os
import json
import time
from sarek_report import make_send_report
from sarek_step import monitor_feedback
from sarek_funcs import compress_results
from sarek_store2db import make_send_result_file

def make_csv_file(analysis_path: str, patient_name: str, sample_type_l: list, sample_name_l: list, sample_file_l: list) -> str:
    csv_file = '{}/samplesheet.csv'.format(analysis_path)
    csv_header = 'patient,sex,status,sample,lane,fastq_1,fastq_2'
    csv_content = ''
    for index in range(0, len(sample_name_l)):
        csv_content += patient_name + ','  # patient
        csv_content += 'XX,'  # sex
        if sample_type_l[index] == 'tumor':  # status
            csv_content += '1,'
        else:
            csv_content += '0,'
        csv_content += sample_name_l[index] + ','  # sample
        csv_content += 'LANE,'  # lane
        csv_content += sample_file_l[index] + '\n'  # fastq1,fastq2
    with open(csv_file, 'w') as csv_f:
        csv_f.write(csv_header + '\n')
        csv_f.write(csv_content[:-1])
    return csv_file


def steward(config_file_path: str, sarek_path: str, genome_base: str, send_message_script: str, sarek_config_d: dict) -> None:
    """
    _summary_

    Args:
        config_file_path (str): _description_
        sarek_path (str): _description_
        genome_base (str): _description_
        send_message_script (str): _description_
    """
    start_time    = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    analysis_path = os.path.dirname(config_file_path)
    with open(config_file_path, 'r') as config_f:
        config_file_d = json.load(config_f)
    # get the data from the config.json
    task_id            = config_file_d['taskId']
    analysis_record_id = config_file_d['analysisRecordId']
    task_name          = config_file_d['taskName']
    pipeline_name      = config_file_d['pipeline']  # panel wes wgs
    patient_name       = config_file_d['patientId2']
    if not patient_name:
        patient_name = '--'
    genome_version     = config_file_d['parameterList']['genome_version']
    mutation_frequency = config_file_d['parameterList']['mutation_site_analysis']
    if mutation_frequency.strip().lower() == 'false':
        is_mutation = False
    else:
        is_mutation = True
    fusion_frequency = config_file_d['parameterList']['gene_fusion_analysis']
    if fusion_frequency.strip().lower() == 'false':
        is_fusion = False
    else:
        is_fusion = True
    if config_file_d['parameterList']['copy_number_variations'].strip().lower() == 'false':
        is_cnv = False
    else:
        is_cnv = True
        if ',' in config_file_d['parameterList']['copy_number_variations']:
            ascat_min_base_qual = config_file_d['parameterList']['copy_number_variations'].split(',')[0]
            ascat_min_counts    = config_file_d['parameterList']['copy_number_variations'].split(',')[1]
            ascat_min_map_qual  = config_file_d['parameterList']['copy_number_variations'].split(',')[2]
            cf_ploidy           = config_file_d['parameterList']['copy_number_variations'].split(',')[3]
        else:
            ascat_min_base_qual = 20
            ascat_min_counts    = 10
            ascat_min_map_qual  = 35
            cf_ploidy = 2
    bed_file = config_file_d['parameterList']['filePath']
    if (bed_file is None) or len(bed_file.strip()) == 0:
        bed_file = sarek_config_d['default_bed_file']
    if 'data_qc' in config_file_d['parameterList']:
        if config_file_d['parameterList']['data_qc'].strip().lower() == 'false':
            data_qc = False
        else:
            data_qc = True
            if ',' in config_file_d['parameterList']['data_qc']:
                split_fastq            = config_file_d['parameterList']['data_qc'].split(',')[0]
                nucleotides_per_second = config_file_d['parameterList']['data_qc'].split(',')[1]
            else:
                split_fastq            = 50000000
                nucleotides_per_second = 1000
    else:
        data_qc                = True
        split_fastq            = 50000000
        nucleotides_per_second = 1000
    if 'align_genome' in config_file_d['parameterList']:
        if config_file_d['parameterList']['align_genome'].strip().lower() == 'false':
            align_genome = False
        else:
            align_genome = True
            aligner      = config_file_d['parameterList']['align_genome'].strip()
    else:
        align_genome = True
        aligner      = 'bwa-mem'
    sample_name_l = []
    sample_type_l = []
    sample_file_l = []
    sample_id_l   = []
    for sample in config_file_d['taskSampleList']:
        sample_name_l.append(sample['sampleName'])  # the sample name
        sample_type_l.append(sample['sampleType'].lower())  # the sample type, [ tumor | normal ]
        sample_file_l.append(sample['fileName'])  # the fastq file including R1 and R2
        sample_id_l.append(sample['sampleId'])
    os.chdir(analysis_path)
    csv_file = make_csv_file(analysis_path, patient_name, sample_type_l, sample_name_l, sample_file_l)
    # make params.json file
    if split_fastq is None:
        split_fastq = 50000000
    if nucleotides_per_second is None:
        nucleotides_per_second = 1000
    if (aligner is None) or (len(aligner) == 0):
        aligner = 'bwa-mem'
    if ascat_min_base_qual is None:
        ascat_min_base_qual = 20
    if ascat_min_counts is None:
        ascat_min_counts = 10
    if ascat_min_map_qual is None:
        ascat_min_map_qual = 35
    if cf_ploidy is None:
        cf_ploidy = 2
    if is_mutation:
        tools = 'strelka,mutect2,freebayes,'
    if is_fusion:
        tools += 'manta,tiddit,'
    if is_cnv:
        tools += 'ascat,cnvkit,controlfreec,'
    if is_mutation or is_fusion or is_cnv:
        tools += 'snpeff'
    params_d = {
        'input'                 : csv_file,
        'outdir'                : 'results',
        'split_fastq'           : split_fastq,
        'nucleotides_per_second': nucleotides_per_second,
        'aligner'               : aligner,
        'igenomes_base'         : genome_base
        }
    if is_mutation or is_fusion or is_cnv:
        params_d['tools'] = tools
    if 'ascat' in tools:
        params_d['ascat_min_base_qual'] = ascat_min_base_qual
        params_d['ascat_min_counts']    = ascat_min_counts
        params_d['ascat_min_map_qual']  = ascat_min_map_qual
    if 'controlfreec' in tools:
        params_d['cf_ploidy'] = cf_ploidy
    if pipeline_name.strip().lower() == 'wes' or pipeline_name.strip().lower() == 'panel':
        params_d['wes']       = True
        params_d['intervals'] = bed_file
    params_file_path = '{}/params.json'.format(analysis_path)
    with open(params_file_path, 'w') as params_f:
        json.dumps(params_d, ensure_ascii=False, fp=params_f, indent=4)
    sarek_command = 'nextflow run -offline -profile singularity -bg -params-file {} {} >>run_sarek.log'.format(params_file_path, sarek_path)
    return_value  = os.system(sarek_command)
    with open('sarek_command.txt', 'w') as sarek_command_f:
        sarek_command_f.write(sarek_command + '\n')
        sarek_command_f.write('return value:{}\n'.format(str(return_value)))
    
    # to monitor the pipeline execution status and send messages
    monitor_feedback_params = {
        'analysis_path'      : analysis_path,
        'send_message_script': send_message_script,
        'call_value'         : return_value,
        'start_time'         : start_time,
        'task_id'            : task_id,
        'analysis_record_id' : analysis_record_id,
        'is_mutation'        : is_mutation,
        'is_fusion'          : is_fusion,
        'is_cnv'             : is_cnv,
        'pipeline_name'      : pipeline_name
    }
    monitor_feedback(monitor_feedback_params)

    # to compress the results for download
    compress_results(analysis_path, '{}/results2download.zip'.format(analysis_path))

    # to generate the report json file
    make_send_report_params = {
        'analysis_path'      : analysis_path,
        'send_message_script': send_message_script,
        'task_id'            : task_id,
        'analysis_record_id' : analysis_record_id,
        'is_mutation'        : is_mutation,
        'is_fusion'          : is_fusion,
        'is_cnv'             : is_cnv,
        'sample_name_l'      : sample_name_l,
        'sample_type_l'      : sample_type_l,
        'pipeline_name'      : pipeline_name
    }
    make_send_report(make_send_report_params)
    
    # to extract the results and store them into the database
    make_send_result_file_params = {
        'analysis_path'      : analysis_path,
        'send_message_script': send_message_script,
        'task_id'            : task_id,
        'analysis_record_id' : analysis_record_id,
        'is_mutation'        : is_mutation,
        'sample_id_l'        : sample_id_l,
        'is_cnv'             : is_cnv,
        'sample_name_l'      : sample_name_l,
        'sample_type_l'      : sample_type_l,
        'genome_version'     : genome_version,
        'pipeline_name'      : pipeline_name
    }
    make_send_result_file(make_send_result_file_params)


def main() -> None:
    parser = argparse.ArgumentParser(description='transfer the config.json to csv file; invoke sarek; feedback the information to front end')
    parser.add_argument('--cfp', required=True, help='the full path for config.json file')
    parser.add_argument('--sarek_path', required=True, help='the full path for sarek')
    parser.add_argument('--genome_base', required=True, help='the full path for genome')
    parser.add_argument('--send_message_script', required=True, help='the full path for the shell script: sendMessage.sh')
    args                = parser.parse_args()
    config_file_path    = args.cfp
    sarek_path          = args.sarek_path
    genome_base         = args.genome_base
    send_message_script = args.send_message_script
    script_path         = os.path.dirname(os.path.abspath(__file__))
    with open('{}/sarek_config.json'.format(script_path), 'r') as sarek_config_f:
        sarek_config_d = json.load(sarek_config_f)
    steward(config_file_path, sarek_path, genome_base, send_message_script, sarek_config_d)

if __name__ == '__main__':
    main()
