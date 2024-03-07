#!/usr/bin/env python3


import argparse
import logging
import os
import json
import time
from sarek_report import make_send_report
from sarek_step import monitor_feedback
from sarek_funcs import compress_results
from sarek_store2db import make_send_result_file

def make_csv_file(csv_file: str, patient_name: str, sample_list: list) -> tuple:
    csv_header       = 'patient,sex,status,sample,lane,fastq_1,fastq_2'
    csv_content      = ''
    all_samples_type = []
    tumor_name       = ''
    normal_name      = ''
    for sample_content in sample_list:
        sex     = 'XX'
        status  = '0' if sample_content['sampleType'] == 'normal' else '1'
        sample  = sample_content['sampleName']
        lane    = 'LANE'
        fastq_1 = sample_content['fileName'].split(',')[0]
        fastq_2 = sample_content['fileName'].split(',')[1]
        csv_content += '{},{},{},{},{},{},{}\n'.format(patient_name, sex, status, sample, lane, fastq_1, fastq_2)
        all_samples_type.append(status)
    with open(csv_file, 'w') as csv_f:
        csv_f.write(csv_header + '\n')
        csv_f.write(csv_content[:-1])
    # to detect the variant call mode
    if ('1' in all_samples_type) and ('0' in all_samples_type):
        variant_calling_mode = 'somatic'
        for sample_content in sample_list:
            if sample_content['sampleType'] == 'tumor':
                tumor_name = sample_content['sampleName']
            elif sample_content['sampleType'] == 'normal':
                normal_name = sample_content['sampleName']
            else:
                pass
    elif ('1' in all_samples_type) and ('0' not in all_samples_type):
        variant_calling_mode = 'tumor_only'
    elif ('1' not in all_samples_type) and ('0' in all_samples_type):
        variant_calling_mode = 'germline'
    else:
        pass
    return variant_calling_mode, tumor_name, normal_name


def steward(config_file_path: str, sarek_path: str, genome_base: str, send_message_script: str, sarek_config_d: dict) -> None:
    """
    _summary_

    Args:
        config_file_path (str): _description_
        sarek_path (str): _description_
        genome_base (str): _description_
        send_message_script (str): _description_
    """
    start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    analysis_path = os.path.dirname(config_file_path)
    os.chdir(analysis_path)
    with open(config_file_path, 'r') as config_f:
        config_file_d = json.load(config_f)
    # get the data from the config.json
    task_id            = config_file_d['taskId']
    analysis_record_id = config_file_d['analysisRecordId']
    task_name          = config_file_d['taskName']
    pipeline_name      = config_file_d['pipeline']  # panel wes wgs
    patient_name       = config_file_d['patientId2'] if config_file_d['patientId2'] else '--'
    genome_version     = config_file_d['parameterList']['genome_version']
    is_mutation        = False if config_file_d['parameterList']['mutation_site_analysis'].strip().lower() == 'false' else True
    is_sv              = False if config_file_d['parameterList']['gene_fusion_analysis'].strip().lower() == 'false' else True
    is_cnv             = False if config_file_d['parameterList']['copy_number_variations'].strip().lower() == 'false' else True
    ## parameters for the ASCAT and Control_FREEC
    if ',' in config_file_d['parameterList']['copy_number_variations']:
        ascat_min_base_qual = config_file_d['parameterList']['copy_number_variations'].split(',')[0]
        ascat_min_counts    = config_file_d['parameterList']['copy_number_variations'].split(',')[1]
        ascat_min_map_qual  = config_file_d['parameterList']['copy_number_variations'].split(',')[2]
        cf_ploidy           = config_file_d['parameterList']['copy_number_variations'].split(',')[3]
    else:
        ascat_min_base_qual = 20
        ascat_min_counts    = 10
        ascat_min_map_qual  = 35
        cf_ploidy           = 2
    bed_file = config_file_d['parameterList']['filePath'] if config_file_d['parameterList']['filePath'] else sarek_config_d['default_bed_file']
    bed_file = bed_file if os.path.isabs(bed_file) else '{}/{}'.format(os.path.dirname(os.path.abspath(__file__)), os.path.basename(bed_file))
    if 'data_qc' in config_file_d['parameterList']:
        data_qc = False if config_file_d['parameterList']['data_qc'].strip().lower() == 'false' else True
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
        align_genome = False if config_file_d['parameterList']['align_genome'].strip().lower() == 'false' else True
        aligner      = config_file_d['parameterList']['align_genome'].strip()
    else:
        align_genome = True
        aligner      = 'bwa-mem'
    # make samplesheet.csv file
    csv_file = '{}/samplesheet.csv'.format(analysis_path)
    variant_calling_mode, tumor_name, normal_name = make_csv_file(csv_file, patient_name, config_file_d['taskSampleList'])
    # make params.json file
    tools = ''
    if variant_calling_mode == 'somatic':
        tools += 'freebayes,strelka,mutect2,' if is_mutation                      else ''
        tools += 'manta,tiddit,'              if is_sv                            else ''
        tools += 'ascat,cnvkit,controlfreec,' if is_cnv                           else ''
        tools += 'snpeff'                     if (is_mutation or is_sv or is_cnv) else ''
    elif variant_calling_mode == 'tumor_only':
        tools += 'freebayes,strelka,mutect2,' if is_mutation                    else ''
        tools += 'manta,tiddit,'              if is_sv                          else ''
        tools += 'cnvkit,controlfreec,'       if is_cnv                         else ''
        tools += 'snpeff'                     if is_mutation or is_sv or is_cnv else ''
    elif variant_calling_mode == 'germline':
        tools += 'freebayes,strelka,' if is_mutation                      else ''
        tools += 'manta,tiddit,'      if is_sv                            else ''
        tools += 'cnvkit,'            if is_cnv                           else ''
        tools += 'snpeff'             if (is_mutation or is_sv or is_cnv) else ''
    else:
        pass
    params_d = {
        'input'                 : csv_file,
        'outdir'                : 'results',
        'split_fastq'           : split_fastq,
        'nucleotides_per_second': nucleotides_per_second,
        'aligner'               : aligner,
        'igenomes_base'         : genome_base
        }
    if is_mutation or is_sv or is_cnv:
        params_d['tools'] = tools
    if 'ascat' in tools:
        params_d['ascat_min_base_qual'] = ascat_min_base_qual
        params_d['ascat_min_counts']    = ascat_min_counts
        params_d['ascat_min_map_qual']  = ascat_min_map_qual
    if 'controlfreec' in tools:
        params_d['cf_ploidy'] = cf_ploidy
    if (pipeline_name.strip().lower() == 'wes') or (pipeline_name.strip().lower() == 'panel'):
        params_d['wes']       = True
        params_d['intervals'] = bed_file
    params_file_path = '{}/params.json'.format(analysis_path)
    with open(params_file_path, 'w') as params_f:
        json.dump(params_d, ensure_ascii=False, fp=params_f, indent=4)
    # to start the sarek
    sarek_command = 'nextflow run -offline -profile singularity -bg -params-file {} {} >>.nextflow_stdout.txt'.format(params_file_path, sarek_path)
    return_value  = os.system(sarek_command)
    logging.info(sarek_command)
    logging.info('return value:{}\n'.format(str(return_value)))
    time.sleep(20)
    # to monitor the pipeline execution status and send messages
    monitor_feedback_params = {
        'analysis_path'      : analysis_path,
        'send_message_script': send_message_script,
        'call_value'         : return_value,
        'start_time'         : start_time,
        'task_id'            : task_id,
        'analysis_record_id' : analysis_record_id,
        'is_mutation'        : is_mutation,
        'is_sv'              : is_sv,
        'is_cnv'             : is_cnv,
        'pipeline_name'      : pipeline_name
    }
    monitor_feedback(monitor_feedback_params)

    # to generate the report json file
    make_send_report_params = {
        'analysis_path'       : analysis_path,
        'send_message_script' : send_message_script,
        'task_id'             : task_id,
        'analysis_record_id'  : analysis_record_id,
        'is_mutation'         : is_mutation,
        'is_sv'               : is_sv,
        'is_cnv'              : is_cnv,
        'tumor_name'          : tumor_name,
        'normal_name'         : normal_name,
        'variant_calling_mode': variant_calling_mode,
        'pipeline_name'       : pipeline_name,
    }
    make_send_report(make_send_report_params)
    
    # to extract the results and store them into the database
    make_send_result_file_params = {
        'analysis_path'       : analysis_path,
        'send_message_script' : send_message_script,
        'task_id'             : task_id,
        'analysis_record_id'  : analysis_record_id,
        'is_mutation'         : is_mutation,
        'is_cnv'              : is_cnv,
        'tumor_name'          : tumor_name,
        'normal_name'         : normal_name,
        'variant_calling_mode': variant_calling_mode,
        'genome_version'      : genome_version,
        'pipeline_name'       : pipeline_name
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
    config_file_path    = config_file_path if os.path.isabs(config_file_path) else os.path.abspath(os.path.basename(config_file_path))
    # logging
    log_file = '{}/sarek.log'.format(os.path.dirname(config_file_path))
    logging.basicConfig(
        filename=log_file, 
        level=logging.INFO, 
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    script_path = os.path.dirname(os.path.abspath(__file__))
    with open('{}/sarek_config.json'.format(script_path), 'r') as sarek_config_f:
        sarek_config_d = json.load(sarek_config_f)
    steward(config_file_path, sarek_path, genome_base, send_message_script, sarek_config_d)

if __name__ == '__main__':
    main()
