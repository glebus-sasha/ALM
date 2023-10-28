#!/bin/bash

# Путь к папке с входными файлами
input_dir="/home/oxkolpakova/data/raw/"

# Путь к папке для сохранения результатов
output_dir="/home/oxkolpakova/data/result/fastp/"

# Параметры для fastp
fastp_params="-q 20 -l 50 --trim_poly_g --thread 12"

# Создаем папку для результатов, если её нет
mkdir -p $output_dir

# Проходим по парам файлов и запускаем fastp для каждой пары
for file1 in ${input_dir}*_R1.fq.gz; do
    file2=${file1/_R1/_R2}  # Получаем имя второго файла (заменяем _R1 на _R2)
    
    # Имя для выходных файлов
    out1="${output_dir}$(basename $file1)"
    out2="${output_dir}$(basename $file2)"
    
    # Имя для HTML-отчета
    report="${output_dir}$(basename $file1 .fq.gz).html"
    
    # Запускаем fastp для текущей пары файлов
    fastp $fastp_params --in1 $file1 --in2 $file2 --out1 $out1 --out2 $out2 -h $report
    
    echo "Обработка файлов $file1 и $file2 завершена. HTML-отчет сохранен в $report"
done

