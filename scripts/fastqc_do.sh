#!/bin/bash

# Путь к папке с входными файлами (результатами fastp)
input_dir="/home/oxkolpakova/data/result/fastp/"

# Путь к папке для сохранения отчетов FastQC
output_dir="/home/oxkolpakova/data/result/fastqc/"

# Создаем папку для отчетов FastQC, если её нет
mkdir -p $output_dir

# Проходим по парам файлов и запускаем FastQC для каждой пары
for file1 in ${input_dir}*_R1.fq.gz; do
    file2=${file1/_R1.fq.gz/_R2.fq.gz}  # Получаем имя второго файла
    
    # Запускаем FastQC для текущей пары файлов
    fastqc -o $output_dir $file1 $file2
    
    echo "FastQC обработал файлы $file1 и $file2. Отчеты сохранены в $output_dir"
done

