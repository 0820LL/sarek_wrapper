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
        if os.path.isdir(os.path.join(input_path, file)):
            get_files(os.path.join(input_path, file), file_list)
        else:
            file_list.append(os.path.join(input_path, file))


def compress_results(intermediate_files_list:list, intermediate_files_zip:str) -> None:
    file_list = []
    for intermediate_files in intermediate_files_list:
        get_files(intermediate_files, file_list)
    with zipfile.ZipFile(intermediate_files_zip, 'w', zipfile.ZIP_DEFLATED) as zip_f:
        for file in file_list:
            zip_f.write(file)
