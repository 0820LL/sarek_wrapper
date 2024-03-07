#!/usr/bin/env python3

import os
import time
from sarek_funcs import send_json_message


def monitor_feedback(params_d:dict) -> None:
    analysis_path       = params_d['analysis_path']
    send_message_script = params_d['send_message_script']
    call_value          = params_d['call_value']
    start_time          = params_d['start_time']
    task_id             = params_d['task_id']
    analysis_record_id  = params_d['analysis_record_id']
    is_mutation         = params_d['is_mutation']
    is_sv               = params_d['is_sv']
    is_cnv              = params_d['is_cnv'],
    pipeline_name       = params_d['pipeline_name']
    os.chdir(analysis_path)
    step_dict = {
        'tTaskId'         : task_id,
        'analysisRecordId': analysis_record_id,
        'pipelineName'    : pipeline_name,
        'analysisStatus'  : '',
        'startDate'       : start_time,
        'endDate'         : '',
        'error'           : 0,
        'taskName'        : 'Step'
    }
    # send the message of step_0
    detect_num = 0
    while True:
        detect_num += 1
        if os.path.exists('{}/results'.format(analysis_path)) or (detect_num > 5):
            break
        else:
            time.sleep(60)
    step_dict['analysisStatus'] = '流程开始'
    step_file_name = 'step_start.json'
    if call_value == 0 and os.path.exists('{}/work'.format(analysis_path)):
        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
        step_dict['error'] = 0
        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
    else:
        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
        step_dict['error'] = 1
        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
        exit('sarek startup failed')
    # send the message of step_x
    while True:
        execution_trace_file = ''
        pipeline_info_dir = '{}/results/pipeline_info'.format(analysis_path)
        if os.path.exists(pipeline_info_dir):
            for file in os.listdir(pipeline_info_dir):
                execution_trace_file = os.path.join(pipeline_info_dir, file) if file.startswith('execution_trace') else ''
                break
        else:
            time.sleep(60)
            continue
        if 'execution_trace' in execution_trace_file:
            break
        else:
            pass
    step_pre = 0
    step_var = 0
    step_fusion = 0
    step_cnv = 0
    step_anno = 0
    step_end = 0
    while step_end == 0:
        # if cancel.txt or Cancel.txt is found, kill the pipeline and exit
        if os.path.exists('{}/cancel.txt'.format(analysis_path)) or os.path.exists('{}/Cancel.txt'.format(analysis_path)):
            step_dict['startDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            step_dict['error'] = 2
            if os.path.exists('{}/.nextflow.pid'.format(analysis_path)):
                with open('{}/.nextflow.pid'.format(analysis_path)) as f:
                    os.system('kill {}'.format(f.read().strip('/n')))
            send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
            exit('The analysis was cancelled')
        with open(execution_trace_file, encoding = 'UTF-8') as trace_f:
            for line in trace_f:
                if step_pre == 0:
                    if 'PREPARE_INTERVALS:CREATE_INTERVALS_BED' in line:
                        step_dict['startDate'] = line.split('\t')[6][:-4]
                        step_dict['analysisStatus'] = '数据预处理'
                        step_file_name = 'step_preprocess.json'
                        continue
                    if ('FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP' in line) and ('COMPLETED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 0
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        step_pre = 1
                        continue
                    if ('FAILED' in line) or ('ABORTED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 1
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        exit('errors in step_preprocess')
                if is_mutation and (step_var == 0):
                    if 'FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP' in line:
                        step_dict['startDate'] = line.split('\t')[6][:-4]
                        step_dict['analysisStatus'] = '突变检测'
                        step_file_name = 'step_callvar.json'
                        continue
                    if ('FILTERMUTECTCALLS' in line) and ('COMPLETED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 0
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        step_var =1
                        continue
                    if ('FAILED' in line) or ('ABORTED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 1
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        exit('errors in step_callvar')
                if is_sv and (step_fusion == 0):
                    if 'TIDDIT_SV' in line:
                        step_dict['startDate'] = line.split('\t')[6][:-4]
                        step_dict['analysisStatus'] = 'SV检测'
                        step_file_name = 'step_sv.json'
                        continue
                    if ('SVDB_MERGE' in line) and ('COMPLETED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 0
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        step_fusion =1
                        continue
                    if ('FAILED' in line) or ('ABORTED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 1
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        exit('errors in step_sv')
                if is_cnv and (step_cnv == 0):
                    if 'CNVKIT_REFERENCE' in line:
                        step_dict['startDate'] = line.split('\t')[6][:-4]
                        step_dict['analysisStatus'] = 'CNV检测'
                        step_file_name = 'step_cnv.json'
                        continue
                    if ('BAM_VARIANT_CALLING_SOMATIC_ALL' in line) and ('CNVKIT_BATCH' in line) and ('COMPLETED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 0
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        step_cnv =1
                        continue
                    if ('FAILED' in line) or ('ABORTED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 1
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        exit('errors in step_cnv')
                if is_mutation and (step_anno == 0):
                    if 'VCF_ANNOTATE_SNPEFF' in line:
                        step_dict['startDate'] = line.split('\t')[6][:-4]
                        step_dict['analysisStatus'] = '突变注释'
                        step_file_name = 'step_anno.json'
                        continue
                    if ('VCF_ANNOTATE_SNPEFF' in line) and ('TABIX_BGZIPTABIX' in line) and ('COMPLETED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 0
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        step_anno = 1
                        continue
                    if ('FAILED' in line) or ('ABORTED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 1
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        exit('errors in step_anno')
            if step_end == 0:
                    if 'MULTIQC' in line:
                        step_dict['startDate'] = line.split('\t')[6][:-4]
                        step_dict['analysisStatus'] = '分析结束'
                        step_file_name = 'step_multiqc.json'
                    if ('MULTIQC' in line) and ('COMPLETED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 0
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        step_end = 1
                        break
                    if ('FAILED' in line) or ('ABORTED' in line):
                        step_dict['endDate'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
                        step_dict['error'] = 1
                        send_json_message(analysis_path, send_message_script, step_dict, step_file_name)
                        exit('errors in step_multiqc')