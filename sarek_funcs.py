#!/usr/bin/env python3

import os
import json
import zipfile

def send_json_message(analysis_path:str, send_message_script: str, message: dict, step_file_name:str) -> None:
    os.system('bash {} \'{}\' 1>>final_message.log 2>&1'.format(send_message_script, json.dumps(message, ensure_ascii=False, indent=4)))
    with open('{}/{}'.format(analysis_path, step_file_name), 'w') as step_f:
        json.dump(message, fp=step_f, ensure_ascii=False, indent=4)


def get_files(input_path:str, file_list:list) -> list:
    files = os.listdir(input_path)
    for file in files:
        if os.path.isdir('{}/{}'.format(input_path, file)):
            get_files('{}/{}'.format(input_path, file), file_list)
        else:
            file_list.append('{}/{}'.format(input_path, file))


def compress_results(analysis_path: str, zip_file:str) -> None:
    file_list = []
    get_files('{}/results/variant_calling'.format(analysis_path), file_list)
    get_files('{}/results/annotation'.format(analysis_path), file_list)
    get_files('{}/results/multiqc'.format(analysis_path), file_list)
    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zip_f:
        for file in file_list:
            zip_f.write(file)
