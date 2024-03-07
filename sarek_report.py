#!/usr/bin/env python3


import os
from sarek_funcs import send_json_message
from sarek_funcs import compress_results


def make_send_report(params_d:dict) -> None:
    analysis_path        = params_d['analysis_path']
    send_message_script  = params_d['send_message_script']
    is_mutation          = params_d['is_mutation']
    is_sv                = params_d['is_sv']
    is_cnv               = params_d['is_cnv']
    variant_calling_mode = params_d['variant_calling_mode']
    tumor_name           = params_d['tumor_name']
    normal_name          = params_d['normal_name']
    pipeline_name        = params_d['pipeline_name']
    
    # to organize the files to download
    multiqc_report = '{}/results/multiqc/multiqc_report.html'.format(analysis_path)
    pipeline_info_dir = '{}/results/pipeline_info'.format(analysis_path)
    execution_report = ''
    for filename in os.listdir(pipeline_info_dir):
        execution_report = os.path.join(pipeline_info_dir, filename) if filename.startswith('execution_report') else ''
    intermediate_files_list = []
    intermediate_files_list.append('{}/results/variant_calling'.format(analysis_path))
    intermediate_files_list.append('{}/results/annotation'.format(analysis_path))
    intermediate_files_list.append('{}/results/reports'.format(analysis_path))
    intermediate_files_zip = '{}/intermediate_files.zip'.format(analysis_path)
    compress_results(intermediate_files_list, intermediate_files_zip)

    report_dict          = {
        'status'          : 'Pass',
        'pipelineName'    : pipeline_name,
        'taskId'          : params_d['task_id'],
        'analysisRecordId': params_d['analysis_record_id'],
        'error'           : 0,
        'taskName'        : 'Report',
        'data'            : []
    }
    report_seq_stat_dict = {
        'sort'     : 1,
        'title'    : '测序数据信息',
        'subtitle1': '测序数据质量评估',
        'memo'     : '',
        'table'    : {
            'data': [
                {
                    'content': '',
                    'path'   : '{}/results/multiqc/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt'.format(analysis_path),
                    'memo'   : '',
                    'preDes' : '使用FastQC软件, 对测序数据进行统计',
                    'title'  : '测序数据统计表',
                    'postDes': ''
                }
            ]
        }, 
        'image': {
            'data': [
                {
                    'title'  : '样本的数据量统计图',
                    'postDes': '其中Duplicate readsd的数量通过软件估计',
                    'memo'   : '',
                    'preDes' : '',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_fastqc_sequence_counts_plot_1.png'.format(analysis_path)
                }
            ]
        },
        'subtitle2': ''
    }
    report_dict['data'].append(report_seq_stat_dict)
    report_seq_qc_dict = {
        'sort'  : 2,
        'title' : '数据质量',
        'memo'  : '',
        'preDes': '数据质控结果展示',
        'image' : {
            'data': [
                {
                    'sort'   : 1,
                    'title'  : '',
                    'postDes': 'reads中每个碱基的平均质量',
                    'memo'   : '',
                    'preDes' : 'Reads的每个碱基平均质量图',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_fastqc_per_base_sequence_quality_plot_1.png'.format(analysis_path)
                }, 
                {
                    'sort'   : 2,
                    'title'  : '',
                    'postDes': '每个Reads的平均质量分布',
                    'memo'   : '',
                    'preDes' : 'Reads的质量值',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_fastqc_per_sequence_quality_scores_plot_1.png'.format(analysis_path)
                }, 
                {
                    'sort'   : 3,
                    'title'  : '',
                    'postDes': '',
                    'memo'   : '每个Reads的平均GC含量',
                    'preDes' : 'Reads的GC含量图',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_fastqc_per_sequence_gc_content_plot_Counts.png'.format(analysis_path)
                }, 
                {
                    'sort'   : 4,
                    'title'  : '',
                    'postDes': '',
                    'memo'   : '',
                    'preDes' : 'FastP过滤reads统计',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_fastp_filtered_reads_plot_1.png'.format(analysis_path)
                }, 
                {
                    'sort'   : 5,
                    'title'  : '',
                    'postDes': '使用GATK picards对duplicate reads进行统计',
                    'memo'   : '',
                    'preDes' : 'Duplicate reads统计',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_picard_deduplication_1.png'.format(analysis_path)
                }
            ]
        }, 
        'subtitle2': ''
    }
    report_dict['data'].append(report_seq_qc_dict)
    report_align_dict = {
        'sort'     : 3,
        'title'    : '比对基因组',
        'subtitle1': '比对结果统计',
        'memo'     : '',
        'table'    : {
            'sort': 1,
            'data': [
                {
                    'content': '',
                    'path'   : '{}/results/multiqc/multiqc_data/mqc_samtools_alignment_plot_1.txt'.format(analysis_path),
                    'memo'   : '',
                    'preDes' : '使用bwa比对参考基因组, samtools stats进行统计',
                    'title'  : '样本比对情况',
                    'postDes': ''
                }
            ]
        }, 
        'image': {
            'data': [
                {
                    'sort'   : 1,
                    'title'  : ' 样本比对统计图',
                    'postDes': '',
                    'memo'   : '',
                    'preDes' : '',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_samtools_alignment_plot_1.png'.format(analysis_path)
                },
                {
                    'sort'   : 2,
                    'title'  : '累计覆盖度统计图',
                    'postDes': '',
                    'memo'   : '',
                    'preDes' : '',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_mosdepth-cumcoverage-dist-id_1.png'.format(analysis_path)
                },
                {
                    'sort'   : 3,
                    'title'  : ' 基因组覆盖度统计图',
                    'postDes': '',
                    'memo'   : '',
                    'preDes' : '',
                    'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_mosdepth-coverage-per-contig_1.png'.format(analysis_path)
                }
            ]
        },
        'subtitle2': ''
    }
    report_dict['data'].append(report_align_dict)
    if is_mutation:
        sort_num = len(report_dict['data']) + 1
        report_mutation_dict = {
            'sort' : sort_num,
            'title': '突变检测',
            'text' : [
                {
                    'content': '',
                    'memo'   : ''
                }
            ],
            'image': {
                'data': [
                    {
                        'sort'   : 1,
                        'title'  : ' 突变替换类型分布图',
                        'postDes': '',
                        'memo'   : '',
                        'preDes' : '使用多种分析工具检测突变',
                        'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_bcftools-stats-subtypes_1.png'.format(analysis_path)
                    },
                    {
                        'sort'   : 2,
                        'title'  : 'Indel长度分布图',
                        'postDes': '',
                        'memo'   : '',
                        'preDes' : '',
                        'path'   : '{}/results/multiqc/multiqc_plots/png/mqc_bcftools_stats_indel-lengths_1.png'.format(analysis_path)
                    }
                ]
            },
            'memo'     : '',
            'subtitle1': '',
            'subtitle2': ''
        }
        report_dict['data'].append(report_mutation_dict)
    if is_cnv:
        sort_num = len(report_dict['data']) + 1
        if variant_calling_mode == 'somatic':
            cnv_file_name = '{}_vs_{}'.format(tumor_name, normal_name)
            cnvkit_result_file = '{0}/results/variant_calling/cnvkit/{1}/{2}-scatter.png'.format(analysis_path, cnv_file_name, tumor_name)
        elif variant_calling_mode == 'germline':
            cnvkit_result_file = '{0}/results/variant_calling/cnvkit/{1}/{2}-scatter.png'.format(analysis_path, normal_name, normal_name)
        elif variant_calling_mode == 'tumor_only':
            cnvkit_result_file = '{0}/results/variant_calling/cnvkit/{1}/{2}-scatter.png'.format(analysis_path, tumor_name, tumor_name)
        else:
            pass
        report_cnv_dict = {
            'sort' : sort_num,
            'title': 'CNV检测',
            'text' : [
                {
                    'content': '',
                    'path'   : '',
                    'memo'   : ''
                }
            ], 
            'image': {
                'data': [
                    {
                        'sort'   : 1,
                        'title'  : ' CNV检测结果展示图',
                        'memo'   : '',
                        'preDes' : '使用CNVkit检测CNV',
                        'path'   : cnvkit_result_file,
                        'postDes': '详细分析结果请下载分析报告'
                    }
                ]
            },
            'memo'     : '',
            'subtitle1': '',
            'subtitle2': ''
        }
        report_dict['data'].append(report_cnv_dict)
    if is_sv:
        report_sv_dict = {
        # there is no picture to display
        }
    sort_num = len(report_dict['data']) + 1
    report_mutation_anno_dict = {
            'sort' : sort_num,
            'title': '突变注释',
            'text' : [
                {
                    'content': '',
                    'path'   : '',
                    'memo'   : ''
                }
            ], 
            'memo' : '',
            'image': {
                'data': []
            },
            'subtitle1': '',
            'subtitle2': ''
        }
    # the file may not exist, need to check
    mqc_snpeff_variant_effects_impact_png = '{}/results/multiqc/multiqc_plots/png/mqc_snpeff_variant_effects_impact_1.png'.format(analysis_path)
    mqc_snpeff_variant_effects_class_png  = '{}/results/multiqc/multiqc_plots/png/mqc_snpeff_variant_effects_class_1.png'.format(analysis_path)
    mqc_snpeff_variant_effects_region_png = '{}/results/multiqc/multiqc_plots/png/mqc_snpeff_variant_effects_region_1.png'.format(analysis_path)
    mqc_snpeff_effects_png                = '{}/results/multiqc/multiqc_plots/png/mqc_snpeff_effects_1.png'.format(analysis_path)
    anno_image_data_sort                  = 1
    if os.path.exists(mqc_snpeff_variant_effects_impact_png):
        mqc_snpeff_variant_effects_impact_png_dict = {
            'sort'   : anno_image_data_sort,
            'title'  : '突变按影响程度分布图',
            'postDes': '',
            'memo'   : '',
            'preDes' : '',
            'path'   : mqc_snpeff_variant_effects_impact_png
        }
        report_mutation_anno_dict['image']['data'].append(mqc_snpeff_variant_effects_impact_png_dict)
        anno_image_data_sort += 1
    if os.path.exists(mqc_snpeff_variant_effects_class_png):
        mqc_snpeff_variant_effects_class_png_dict = {
            'sort'   : anno_image_data_sort,
            'title'  : '突变按功能分布图',
            'postDes': '',
            'memo'   : '',
            'preDes' : '',
            'path'   : mqc_snpeff_variant_effects_class_png
        }
        report_mutation_anno_dict['image']['data'].append(mqc_snpeff_variant_effects_class_png_dict)
        anno_image_data_sort += 1
    if os.path.exists(mqc_snpeff_variant_effects_region_png):
        mqc_snpeff_variant_effects_region_png_dict = {
            'sort'   : anno_image_data_sort,
            'title'  : '突变按基因组区域分布图',
            'postDes': '',
            'memo'   : '',
            'preDes' : '',
            'path'   : mqc_snpeff_variant_effects_region_png
        }
        report_mutation_anno_dict['image']['data'].append(mqc_snpeff_variant_effects_region_png_dict)
        anno_image_data_sort += 1
    if os.path.exists(mqc_snpeff_effects_png):
        mqc_snpeff_effects_png_dict = {
            'sort'   : anno_image_data_sort,
            'title'  : '突变影响类型分布图',
            'postDes': '',
            'memo'   : '',
            'preDes' : '',
            'path'   : mqc_snpeff_effects_png
        }
        report_mutation_anno_dict['image']['data'].append(mqc_snpeff_effects_png_dict)
    report_dict['data'].append(report_mutation_anno_dict)
    sort_num = len(report_dict['data']) + 1
    report_reference_dict = {
        'sort' : sort_num,
        'title': '参考文献',
        'text' : [
            {
                "sort"   : 1,
                "content": "[1] Garcia MU, Juhos S, Larsson M, Olason PI, Martin M, Eisfeldt J, DiLorenzo S, Sandgren J, Díaz De Ståhl T, Ewels PA, Wirta V, Nistér M, Käller M, Nystedt B. Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants. F1000Res. 2020 Jan 29;9:63. eCollection 2020. doi: 10.12688/f1000research.16665.2. PubMed PMID: 32269765.",
                "memo"   : ""
            },
            {
                "sort"   : 2,
                "content": "[2] Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.",
                "memo"   : ""
            },
            {
                "sort"   : 3,
                "content": "[3] Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.",
                "memo"   : ""
            },
            {
                "sort"   : 4,
                "content": "[4] McKenna A, Hanna M, Banks E, et al.: The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. doi: 10.1101/gr.107524.110. Epub 2010 Jul 19. PubMed PMID: 20644199; PubMed Central PMCID: PMC2928508.",
                "memo"   : ""
            },
            {
                "sort"   : 5,
                "content": "[5] Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PubMed PMID: 27312411. PubMed Central PMCID: PMC5039924.",
                "memo"   : ""
            },
            {
                "sort"   : 6,
                "content": "[6] Cingolani P, Platts A, Wang le L, et al.: A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly (Austin). Apr-Jun 2012;6(2):80-92. doi: 10.4161/fly.19695. PubMed PMID: 22728672; PubMed Central PMCID: PMC3679285.",
                "memo"   : ""
            },
            {
                "sort"   : 7,
                "content": "[7] Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.",
                "memo"   : ""
            },
            {
                "sort"   : 8,
                "content": "[8] Li H: Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv 2013. doi: 10.48550/arXiv.1303.3997",
                "memo"   : ""
            },
            {
                "sort"   : 9,
                "content": "[9] Li H: A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. doi: 10.1093/bioinformatics/btr509. PubMed PMID: 21903627; PubMed Central PMCID: PMC3198575.",
                "memo"   : ""
            }
        ]
    }
    report_dict['data'].append(report_reference_dict)
    sort_num = len(report_dict['data']) + 1
    report_download_dict = {
        'sort' : sort_num,
        'title': '结果文件下载',
        'text' : [
            {
                'sort'   : 1,
                'title'  : 'MultiQC报告下载',
                'content': 'MultiQC下载：#&{}'.format(multiqc_report),
                'postDes': '详细的分析结果请下载并查看',
                'memo'   : '',
                'preDes' : ''
            },
            {
                'sort'   : 2,
                'title'  : 'pipeline运行报告下载',
                'content': 'pipelines运行报告下载：#&{}'.format(execution_report),
                'postDes': '详细的分析分析结果请下载并查看',
                'memo'   : '',
                'preDes' : ''
            },
            {
                'sort'   : 3,
                'title'  : '分析过程中间文件下载',
                'content': '分析过程中间文件下载：#&{}'.format(intermediate_files_zip),
                'postDes': '详细的分析结果请下载查看',
                'memo'   : '',
                'preDes' : ''
            }
        ], 
        'memo'     : '',
        'subtitle1': '分析报告结果下载',
        'subtitle2': ''
    }
    report_dict['data'].append(report_download_dict)
    send_json_message(analysis_path, send_message_script, report_dict, 'Report.json')