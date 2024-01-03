#!/usr/bin/env python3 


import argparse
import gzip
import json
import os
from sarek_funcs import send_json_message


def get_specified_suffix_file(files: list, suffix: str) -> str:
    for file in files:
        if file.lower().endswith(suffix):
            return file


def write_qc_table(task_id: int, analysis_record_id: int, sample_name_l: list, analysis_path: str) -> str:
    '''
    _summary_

    Args:
        task_id (int)           : _description_
        analysis_record_id (int): _description_
        sample_name_l (list)    : _description_
        analysis_path (str)     : _description_
    '''
    qc_table_d = {}
    for sample in sample_name_l:
        qc_table_d[sample] = {}
        # extract some qc data from the fastp results
        fastp_result_path = '{}/results/reports/fastp/{}'.format(
            analysis_path, sample)
        fastp_json_file = fastp_result_path + '/' + \
            get_specified_suffix_file(os.listdir(
                fastp_result_path), '.fastp.json')
        with open(fastp_json_file, 'r') as fastp_json_f:
            fastp_d = json.load(fastp_json_f)
        qc_table_d[sample]['total_reads_before']       = fastp_d['summary']['before_filtering']['total_reads']
        qc_table_d[sample]['total_bases_before']       = fastp_d['summary']['before_filtering']['total_bases']
        qc_table_d[sample]['q20_bases_before']         = fastp_d['summary']['before_filtering']['q20_bases']
        qc_table_d[sample]['q30_bases_before']         = fastp_d['summary']['before_filtering']['q30_bases']
        qc_table_d[sample]['q20_rate_before']          = fastp_d['summary']['before_filtering']['q20_rate']
        qc_table_d[sample]['q30_rate_before']          = fastp_d['summary']['before_filtering']['q30_rate']
        qc_table_d[sample]['read1_mean_length_before'] = fastp_d['summary']['before_filtering']['read1_mean_length']
        qc_table_d[sample]['read2_mean_length_before'] = fastp_d['summary']['before_filtering']['read2_mean_length']
        qc_table_d[sample]['gc_content_before']        = fastp_d['summary']['before_filtering']['gc_content']
        qc_table_d[sample]['total_reads_after']        = fastp_d['summary']['after_filtering']['total_reads']
        qc_table_d[sample]['total_bases_after']        = fastp_d['summary']['after_filtering']['total_bases']
        qc_table_d[sample]['q20_bases_after']          = fastp_d['summary']['after_filtering']['q20_bases']
        qc_table_d[sample]['q30_bases_after']          = fastp_d['summary']['after_filtering']['q30_bases']
        qc_table_d[sample]['q20_rate_after']           = fastp_d['summary']['after_filtering']['q20_rate']
        qc_table_d[sample]['q30_rate_after']           = fastp_d['summary']['after_filtering']['q30_rate']
        qc_table_d[sample]['read1_mean_length_after']  = fastp_d['summary']['after_filtering']['read1_mean_length']
        qc_table_d[sample]['read2_mean_length_after']  = fastp_d['summary']['after_filtering']['read2_mean_length']
        qc_table_d[sample]['gc_content_after']         = fastp_d['summary']['after_filtering']['gc_content']
        qc_table_d[sample]['duplication_rate']         = fastp_d['duplication']['rate']
        # extract some map related data from the samtools stats results
        samtools_stats_result_file = '{0}/results/reports/samtools/{1}/{1}.md.cram.stats'.format(
            analysis_path, sample)
        samtools_stats_d = {}
        with open(samtools_stats_result_file, 'r') as samtools_stats_f:
            for line in samtools_stats_f:
                if line.startswith('SN'):
                    line_l                                 = line[:-1].split('\t')
                    samtools_stats_d[line_l[1].strip(':')] = line_l[2]
        qc_table_d[sample]['mapped_bases']       = samtools_stats_d['bases mapped (cigar)']
        qc_table_d[sample]['unmapped_bases']     = str(int(qc_table_d[sample]['total_bases_after']) - int(qc_table_d[sample]['mapped_bases']))
        qc_table_d[sample]['mapped_reads']       = samtools_stats_d['reads mapped']
        qc_table_d[sample]['unmapped_reads']     = samtools_stats_d['reads unmapped']
        qc_table_d[sample]['duplicated_reads']   = samtools_stats_d['reads duplicated']
        qc_table_d[sample]['alignment_rate']     = str(round(int(samtools_stats_d['reads mapped']) / int(samtools_stats_d['sequences']), 4))
        qc_table_d[sample]['mean_read_length']   = samtools_stats_d['average length']
        qc_table_d[sample]['insert_size_median'] = samtools_stats_d['insert size average']
        qc_table_d[sample]['error_rate']         = samtools_stats_d['error rate']
        # extract some coverage data from the mosdepth results
        mosdepth_result_file = '{0}/results/reports/mosdepth/{1}/{1}.recal.mosdepth.summary.txt'.format(
            analysis_path, sample)
        with open(mosdepth_result_file, 'r') as mosdepth_result_f:
            for line in mosdepth_result_f:
                if line.startswith('total_region'):
                    qc_table_d[sample]['target_mean_coverage'] = line.split('\t')[3]
        # set as 'NA' for the metrics that can not find in the sarek results
        qc_table_d[sample]['file_name']                  = 'NA'
        qc_table_d[sample]['on_target_reads']            = 'NA'
        qc_table_d[sample]['off_target_reads']           = 'NA'
        qc_table_d[sample]['capture_efficiency']         = 'NA'
        qc_table_d[sample]['fraction_of_target_regions'] = 'NA'
        qc_table_d[sample]['uniformity_of_coverage']     = 'NA'
        qc_table_d[sample]['mean_mapping_quality']       = 'NA'
        qc_table_d[sample]['homopolyer_indels']          = 'NA'
    # prepare the table content and write the QC metrics into a file
    table_header  = 'taskId\t'
    table_header += 'analysisRecordId\t'
    table_header += 'sampleName\t'
    table_header += 'fileName\t'
    table_header += 'totalBases\t'
    table_header += 'mappedBases\t'
    table_header += 'unmappedBases\t'
    table_header += 'q20\t'
    table_header += 'q30\t'
    table_header += 'totalReads\t'
    table_header += 'mappedReads\t'
    table_header += 'unmappedReads\t'
    table_header += 'alignmentRate\t'
    table_header += 'meanReadLength\t'
    table_header += 'onTargetReads\t'
    table_header += 'offTargetReads\t'
    table_header += 'captureEfficiency\t'
    table_header += 'duplicationRate\t'
    table_header += 'targetMeanCoverage\t'
    table_header += 'fractionOfTargetRegions\t'
    table_header += 'gcPercentage\t'
    table_header += 'uniformityOfCoverage\t'
    table_header += 'meanMappingQuality\t'
    table_header += 'insertSizeMedian\t'
    table_header += 'errorRate\t'
    table_header += 'homopolyerIndels'
    table_content = ''
    for sample in sample_name_l:
        table_content += str(task_id) + '\t'
        table_content += str(analysis_record_id) + '\t'
        table_content += sample + '\t'  # sampleName
        table_content += qc_table_d[sample]['file_name'] + '\t'  # fileName
        # totalBases
        table_content += str(qc_table_d[sample]['total_bases_after']) + '\t'
        # mappedBases
        table_content += qc_table_d[sample]['mapped_bases'] + '\t'
        # unmappedBases
        table_content += qc_table_d[sample]['unmapped_bases'] + '\t'
        table_content += str(qc_table_d[sample]
                            ['q20_rate_after']) + '\t'  # q20
        table_content += str(qc_table_d[sample]
                            ['q30_rate_after']) + '\t'  # q30
        # totalReads
        table_content += str(qc_table_d[sample]['total_reads_after']) + '\t'
        # mappedReads
        table_content += qc_table_d[sample]['mapped_reads'] + '\t'
        # unmappedReads
        table_content += qc_table_d[sample]['unmapped_reads'] + '\t'
        # alignmentRate
        table_content += qc_table_d[sample]['alignment_rate'] + '\t'
        # meanReadLength
        table_content += qc_table_d[sample]['mean_read_length'] + '\t'
        # onTargetReads
        table_content += qc_table_d[sample]['on_target_reads'] + '\t'
        # offTargetReads
        table_content += qc_table_d[sample]['off_target_reads'] + '\t'
        # captureEfficiency
        table_content += qc_table_d[sample]['capture_efficiency'] + '\t'
        # duplicationRate
        table_content += str(qc_table_d[sample]['duplication_rate']) + '\t'
        # targetMeanCoverage
        table_content += qc_table_d[sample]['target_mean_coverage'] + '\t'
        # fractionOfTargetRegions
        table_content += qc_table_d[sample]['fraction_of_target_regions'] + '\t'
        # gcPercentage
        table_content += str(qc_table_d[sample]['gc_content_after']) + '\t'
        # uniformityOfCoverage
        table_content += qc_table_d[sample]['uniformity_of_coverage'] + '\t'
        # meanMappingQuality
        table_content += qc_table_d[sample]['mean_mapping_quality'] + '\t'
        # insertSizeMedian
        table_content += qc_table_d[sample]['insert_size_median'] + '\t'
        table_content += qc_table_d[sample]['error_rate'] + '\t'  # errorRate
        # homopolyerIndels
        table_content += qc_table_d[sample]['homopolyer_indels'] + '\n'
    qc_table_file = '{}/table_qc.txt'.format(analysis_path)
    with open(qc_table_file, 'w') as qc_table_f:
        qc_table_f.write(table_header + '\n')
        qc_table_f.write(table_content[:-1])
    return qc_table_file


def write_variants_table(task_id: int, analysis_record_id: int, sample_name_l: list, sample_id_l: list, sample_type_l: list, genome_version: str, analysis_path: str) -> str:
    '''
    extract variants data from the mutect2 results: tumor_vs_normal

    Args:
        task_id (int)           : _description_
        analysis_record_id (int): _description_
        sample_name_l (list)    : _description_
        sample_id_l (list)      : _description_
        sample_type_l (list)    : _description_
        genome_version (str)    : _description_
        analysis_path (str)     : _description_
    '''
    tumor_index           = sample_type_l.index('tumor')
    normal_index          = sample_type_l.index('normal')
    tumor_vs_normal       = '{}_vs_{}'.format(sample_name_l[tumor_index], sample_name_l[normal_index])
    tumor_id_vs_normal_id = '{}_vs_{}'.format(sample_id_l[tumor_index], sample_id_l[normal_index])
    mutect2_result_file   = '{0}/results/annotation/mutect2/{1}/{1}.mutect2.filtered_snpEff.ann.vcf.gz'.format(analysis_path, tumor_vs_normal)
    mutect2_result_l      = []
    with gzip.open(mutect2_result_file, 'rt') as mutect2_result_f:
        for line in mutect2_result_f:
            if line.startswith('chr') and line.split('\t')[6] == 'PASS':
                mutect2_result_l.append(line[:-1])
    table_header   = 'taskId\t'
    table_header  += 'analysisRecordId\t'
    table_header  += 'sampleName\t'
    table_header  += 'sampleId\t'
    table_header  += 'chromosome\t'
    table_header  += 'hgPosition\t'
    table_header  += 'ref\t'
    table_header  += 'alt\t'
    table_header  += 'aminoAcidChange\t'
    table_header  += 'geneName\t'
    table_header  += 'geneBankId\t'
    table_header  += 'aminoAcidLength\t'
    table_header  += 'pValue\t'
    table_header  += 'codonChange\t'
    table_header  += 'effect\t'
    table_header  += 'exonRank\t'
    table_header  += 'functionalClass\t'
    table_header  += 'geneCoding\t'
    table_header  += 'transcriptId\t'
    table_header  += 'cosmic\t'
    table_header  += 'resultsFilePath\t'
    table_header  += 'varfraqPrecentMax\t'
    table_header  += 'covMax\t'
    table_header  += 'hgVerson\t'
    table_header  += 'ci\t'
    table_header  += 'genomeType\t'
    table_header  += 'tag\t'
    table_header  += 'snpId\t'
    table_header  += 'polyphen2HvarScore\t'
    table_header  += 'siftScore\t'
    table_header  += 'globalMaf\t'
    table_header  += 'raceMaf\t'
    table_header  += 'altPair\t'
    table_header  += 'lociPair\t'
    table_header  += 'Phase\t'
    table_header  += 'phaseFreq\t'
    table_header  += 'aberrationId\t'
    table_header  += 'source\t'
    table_header  += 'cancerType\t'
    table_header  += 'geneticAlteration\t'
    table_header  += 'diseaseId'
    table_content  = ''
    for mutect2_result in mutect2_result_l:
        mutect2_result_split  = mutect2_result.split('\t')
        mutect2_info_split    = mutect2_result_split[7].split('|')
        table_content        += str(task_id) + '\t'  # taskId
        table_content        += str(analysis_record_id) + '\t'  # analysisRecordId
        table_content        += tumor_vs_normal + '\t'  # sampleName
        table_content        += tumor_id_vs_normal_id + '\t'  # sampleId
        if mutect2_info_split[0]:
            table_content += mutect2_result_split[0] + '\t'  # chromosome
        else:
            table_content += 'NA\t'  # chromosome
        if mutect2_result_split[1]:
            table_content += mutect2_result_split[1] + '\t'  # hgPosition
        else:
            table_content += 'NA\t'  # hgPosition
        if mutect2_result_split[3]:
            table_content += mutect2_result_split[3] + '\t'  # ref
        else:
            table_content += 'NA\t'  # ref
        if mutect2_result_split[4]:
            table_content += mutect2_result_split[4] + '\t'  # alt
        else:
            table_content += 'NA\t'  # alt
        if mutect2_info_split[11]:
            table_content += mutect2_info_split[11] + '\t'  # aminoAcidChange
        else:
            table_content += 'NA\t'  # aminoAcidChange
        if mutect2_info_split[4]:
            table_content += mutect2_info_split[4] + '\t'  # geneName
        else:
            table_content += 'NA\t'  # geneName
        if mutect2_info_split[5]:
            table_content += mutect2_info_split[5] + '\t'  # geneBankId
        else:
            table_content += 'NA\t'  # geneBankId
        table_content += 'NA\t'  # aminoAcidLength
        table_content += 'NA\t'  # pValue
        table_content += 'NA\t'  # codonChange
        if mutect2_info_split[2]:
            table_content += mutect2_info_split[2] + '\t'  # effect
        else:
            table_content += 'NA\t'  # effect
        table_content += 'NA\t'  # exonRank
        if mutect2_info_split[3]:
            table_content += mutect2_info_split[3] + '\t'  # functionalClass
        else:
            table_content += 'NA\t'  # functionalClass
        table_content += 'NA\t'  # geneCoding
        if mutect2_info_split[7]:
            table_content += mutect2_info_split[7] + '\t'  # transcriptId
        else:
            table_content += 'NA\t'  # transcriptId
        table_content += 'NA\t'  # cosmic
        table_content += mutect2_result_file + '\t'  # resultsFilePath
        table_content += 'NA\t'  # varfraqPrecentMax
        table_content += 'NA\t'  # covMax
        table_content += genome_version + '\t'  # hgVerson
        table_content += 'NA\t'  # ci
        table_content += 'NA\t'  # genomeType
        table_content += 'NA\t'  # tag
        table_content += 'NA\t'  # snpId
        table_content += 'NA\t'  # polyphen2HvarScore
        table_content += 'NA\t'  # siftScore
        table_content += 'NA\t'  # globalMaf
        table_content += 'NA\t'  # raceMaf
        table_content += 'NA\t'  # altPair
        table_content += 'NA\t'  # lociPair
        table_content += 'NA\t'  # Phase
        table_content += 'NA\t'  # phaseFreq
        table_content += 'NA\t'  # aberrationId
        table_content += 'NA\t'  # source
        table_content += 'NA\t'  # cancerType
        table_content += 'NA\t'  # geneticAlteration
        table_content += 'NA\n'  # diseaseId
    variants_table_file = '{}/table_variants.txt'.format(analysis_path)
    with open(variants_table_file, 'w') as variants_table_f:
        variants_table_f.write(table_header + '\n')
        variants_table_f.write(table_content[:-1])
    return variants_table_file


def write_cnv_table(task_id: int, analysis_record_id: int, sample_name_l: list, sample_id_l: list, sample_type_l: list, genome_version: str, analysis_path: str) -> str:
    '''
    extract cnv data from the cnvkit results: tumor_vs_normal

    Args:
        task_id (int)           : _description_
        analysis_record_id (int): _description_
        sample_name_l (list)    : _description_
        sample_id_l (list)      : _description_
        sample_type_l (list)    : _description_
        genome_version (str)    : _description_
        analysis_path (str)     : _description_
    '''
    tumor_index = sample_type_l.index('tumor')
    normal_index = sample_type_l.index('normal')
    tumor_vs_normal = '{}_vs_{}'.format(
        sample_name_l[tumor_index], sample_name_l[normal_index])
    tumor_id_vs_normal_id = '{}_vs_{}'.format(
        sample_id_l[tumor_index], sample_id_l[normal_index])
    cnvkit_result_file = '{0}/results/variant_calling/cnvkit/{1}/{2}.call.cns'.format(
        analysis_path, tumor_vs_normal, sample_name_l[tumor_index])
    cnvkit_result_l = []
    with open(cnvkit_result_file, 'rt') as cnvkit_result_f:
        for line in cnvkit_result_f:
            if line.startswith('chr'):
                cnvkit_result_l.append(line[:-1])
    table_header = 'taskId\t'
    table_header += 'analysisRecordId\t'
    table_header += 'sampleName\t'
    table_header += 'sampleId\t'
    table_header += 'chromosome\t'
    table_header += 'copyNumber\t'
    table_header += 'geneName\t'
    table_header += 'hgPosiStart\t'
    table_header += 'hgPosiEnd\t'
    table_header += 'hgVerson\t'
    table_header += 'aberrationId\t'
    table_header += 'source\t'
    table_header += 'cancerType\t'
    table_header += 'geneticAlteration\t'
    table_header += 'diseaseId'
    table_content = ''
    for cnvkit_result in cnvkit_result_l:
        if cnvkit_result.startswith('chromosome'):
            continue
        cnvkit_result_split  = cnvkit_result.split('\t')
        table_content       += str(task_id) + '\t'  # taskId
        table_content       += str(analysis_record_id) + '\t'  # analysisRecordId
        table_content       += tumor_vs_normal + '\t'  # sampleName
        table_content       += tumor_id_vs_normal_id + '\t'  # sampleId
        table_content       += cnvkit_result_split[0] + '\t'  # chromosome
        table_content       += cnvkit_result_split[5] + '\t'  # copyNumber
        gene_name            = cnvkit_result_split[3]
        if gene_name.startswith('-'):
            gene_name = gene_name.split(',')[1].split('|')[1]
        else:
            gene_name = gene_name.split(',')[0].split('|')[1]
        table_content += gene_name + '\t'  # geneName
        table_content += cnvkit_result_split[1] + '\t'  # hgPosiStart
        table_content += cnvkit_result_split[2] + '\t'  # hgPosiEnd
        table_content += genome_version + '\t'  # hgVerson
        table_content += 'NA\t'  # aberrationId
        table_content += 'NA\t'  # source
        table_content += 'NA\t'  # cancerType
        table_content += 'NA\t'  # geneticAlteration
        table_content += 'NA\n'  # diseaseId
    cnv_table_file = '{}/table_cnv.txt'.format(analysis_path)
    with open(cnv_table_file, 'w') as cnv_table_f:
        cnv_table_f.write(table_header + '\n')
        cnv_table_f.write(table_content[:-1])
    return cnv_table_file


def extract_fusion_data():
    pass

def main2(analysis_path: str, config_file_path: str) -> None:
    with open(config_file_path, 'r') as config_f:
        config_file_d = json.load(config_f)
    # get the data in the config.json
    task_id            = config_file_d['taskId']
    analysis_record_id = config_file_d['analysisRecordId']
    task_name          = config_file_d['taskName']
    pipeline_name      = config_file_d['pipeline']
    patient_name       = config_file_d['patientId2']
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
    if config_file_d['parameterList']['copy_number_variations'].strip().lower() == 'true': # [ true | false ]
        is_cnv = True
    else:
        is_cnv = False
    # get the sample names from the file config.json and then store it in a list
    sample_name_l = []
    sample_id_l   = []
    sample_type_l = []
    for task_sample in config_file_d['taskSampleList']:
        sample_name_l.append(task_sample['sampleName'])
        sample_id_l.append(task_sample['sampleId'])
        sample_type_l.append(task_sample['sampleType'].lower())
    # extract the QC metrics and write them to a file
    qc_table_file = write_qc_table(task_id, analysis_record_id, sample_name_l, analysis_path)
    # extract the variants information from the mutect2 vcf file and write them to a file
    if is_mutation:
        variants_table_file = write_variants_table(task_id, analysis_record_id, sample_name_l,
                            sample_id_l, sample_type_l, genome_version, analysis_path)
    # extract the cnv information from the cnvkit result file and write them to a file
    if is_cnv:
        cnv_table_file = write_cnv_table(task_id, analysis_record_id, sample_name_l,
                        sample_id_l, sample_type_l, genome_version, analysis_path)


def make_send_result_file(params_d: dict) -> None:
    task_id             = params_d['task_id']
    analysis_record_id  = params_d['analysis_record_id']
    analysis_path       = params_d['analysis_path']
    is_mutation         = params_d['is_mutation']
    is_cnv              = params_d['is_cnv']
    sample_id_l         = params_d['sample_id_l']
    sample_name_l       = params_d['sample_name_l']
    sample_type_l       = params_d['sample_type_l']
    genome_version      = params_d['genome_version']
    send_message_script = params_d['send_message_script']
    pipeline_name       = params_d['pipeline_name']
    result_dict = {
        'status'          : 'Pass',
        'pipelineName'    : pipeline_name,
        'taskId'          : task_id,
        'analysisRecordId': analysis_record_id,
        'error'           : 0,
        'taskName'        : 'Result',
        'qc_table'        : None,
        'variant_table'   : None,
        'cnv_table'       : None,
        'fusion_table'    : None,
        'as_table'        : None,
        'de_table'        : None
    }
    qc_table_file           = write_qc_table(task_id, analysis_record_id, sample_name_l, analysis_path)
    result_dict['qc_table'] = qc_table_file
    if is_mutation:
        variants_table_file          = write_variants_table(task_id, analysis_record_id, sample_name_l, sample_id_l, sample_type_l, genome_version, analysis_path)
        result_dict['variant_table'] = variants_table_file
    if is_cnv:
        cnv_table_file           = write_cnv_table(task_id, analysis_record_id, sample_name_l, sample_id_l, sample_type_l, genome_version, analysis_path)
        result_dict['cnv_table'] = cnv_table_file
    send_json_message(analysis_path, send_message_script, result_dict, 'Results.json')


def main() -> None:
    parser = argparse.ArgumentParser(description='extract data and transfer format for the mysql database')
    parser.add_argument('--analysis_path', required=True, help='the results path of certain analysis')
    parser.add_argument('--cfp', required=True, help='the full path for config.json file')
    args = parser.parse_args()
    analysis_path = args.analysis_path
    config_file_path = args.cfp
    main2(analysis_path, config_file_path)


if __name__ == '__main__':
    main()