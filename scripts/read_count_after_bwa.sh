# Проверяем существование папки /home/oxkolpakova/summary и создаем ее, если не существует
summary_dir="/home/oxkolpakova/summary"
if [ ! -d "$summary_dir" ]; then
    mkdir -p "$summary_dir"
fi

# Создаем файл для записи результатов
output_file="${summary_dir}/num_of_reads_after_bwa"

# Очищаем или создаем файл с результатами
> "$output_file"

# Создаем заголовок таблицы
echo -e "номер файла\tчисло ридов в sorted\tчисло ридов в unmapped\tсумма ридов\tпроцент unmapped от суммы" > "$output_file"

# Инициализируем переменные для суммы ридов и процента
total_reads=0

# Запускаем цикл для файлов с 10 по 22
for i in {1..22}; do
    # Формируем путь к файлам
    sorted_file="/home/oxkolpakova/data/result/bwa/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_${i}_L00_sorted.bam"
    unmapped_file="/home/oxkolpakova/data/result/bwa/unmapped/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_${i}_L00.bam"
    
    # Подсчитываем число ридов в каждом файле
    sorted_count=$(samtools view -c "$sorted_file")
    unmapped_count=$(samtools view -c "$unmapped_file")
    
    # Вычисляем сумму ридов
    total_reads=$((total_reads + sorted_count + unmapped_count))
    
    # Проверяем, что total_reads больше нуля перед вычислением процента
    if [ "$total_reads" -gt 0 ]; then
        percent_unmapped=$(bc <<< "scale=2; 100 * $unmapped_count / $total_reads")
    else
        percent_unmapped="N/A" # Если total_reads равно нулю, ставим "N/A" для процента
    fi

    # Выводим результат в виде строки таблицы
    echo -e "${i}\t${sorted_count}\t${unmapped_count}\t${total_reads}\t${percent_unmapped}%" >> "$output_file"
done

