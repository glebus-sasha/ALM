#!/bin/bash

# Полный путь к .sif файлу
sif_path='/home/oxkolpakova/programs/gapseq_latest.sif'

# Директория для вывода логов
log_dir='/home/oxkolpakova/ALM/scripts/logs/'

# Путь для файла media.csv
media_csv_path='/home/oxkolpakova/ALM/data/MM_anaerobic_Acetate_H2.csv'

# Функция для обработки одного файла
process_file() {
    local input_file=$1
    
    # Имя файла без расширения
    xbase=$(basename $input_file)
    shname=${xbase%.*}
    
    # Журнал для вывода
    log_file="$log_dir/$shname.log"

    # Запустить последние две gapseq команды внутри контейнера и логгировать
    #singularity run -B $output_dir:/mnt $sif_path gapseq find -p all -b 100 -m Bacteria /mnt/$xbase > $log_file 2>&1
    #singularity run -B $output_dir:/mnt $sif_path gapseq find-transport -b 100 /mnt/$xbase >> $log_file 2>&1
    #singularity run -B $log_dir:/mnt $sif_path gapseq draft -r /home/oxkolpakova/ALM/scripts/$shname-all-Reactions.tbl -t /home/oxkolpakova/ALM/scripts/$shname-Transporter.tbl -p /home/oxkolpakova/ALM/scripts/$shname-all-Pathways.tbl -c $input_file -u 100 -l 50 > $log_file 2>&1
    singularity run -B $log_dir:/mnt $sif_path gapseq fill -m /home/oxkolpakova/ALM/scripts/$shname-draft.RDS -n $media_csv_path -c /home/oxkolpakova/ALM/scripts/$shname-rxnWeights.RDS -b 100 -g /home/oxkolpakova/ALM/scripts/$shname-rxnXgenes.RDS >> $log_file 2>&1

    echo "Обработка файла $input_file завершена. Журнал: $log_file"
}

# Массив файлов для обработки
files=('/home/oxkolpakova/ALM/data/results/meta_SPAdes/metaSPAdes_S1_Contigs.fasta'
       '/home/oxkolpakova/ALM/data/results/meta_SPAdes/metaSPAdes_S2_Contigs.fasta'
       '/home/oxkolpakova/ALM/data/results/meta_SPAdes/metaSPAdes_S3_Contigs.fasta'
       '/home/oxkolpakova/ALM/data/results/meta_SPAdes/metaSPAdes_S4_Contigs.fasta')

# Создать директорию для логов, если она не существует
mkdir -p $log_dir

# Запустить обработку каждого файла в фоне
for file in "${files[@]}"; do
    process_file "$file" &
done

# Ожидать завершения всех фоновых процессов
wait

echo "Все файлы обработаны."

