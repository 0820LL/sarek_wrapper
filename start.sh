#!/usr/bin/bash

if [ "$#" -ne 1 ]; then
    echo "The script need exactly 1 argument: config.json"
    echo "/.../$0 /.../config.json"
    exit
fi

config_file=$(realpath "$1")
analysis_dir=$(dirname "$config_file")
script_path=$(dirname "$(realpath "$0")")
configuration_file="$script_path"/../../config/configuration.json
genome_base=$(jq ".genome_base" "$configuration_file" | sed 's/\"//g')
sendMessage=$(jq ".jms" "$configuration_file" | sed 's/\"//g')

cd "$analysis_dir" || exit
python3 "$script_path"/sarek_start.py \
    --cfp "$config_file" \
    --sarek_path "$script_path"/../nf-core-sarek-3.1.2/workflow/main.nf \
    --genome_base "$genome_base" \
    --send_message_script "$sendMessage"
