#!/bin/bash

NUM_OF_RUNS=100

if [ ! -z "$1" ]; then
    NUM_OF_RUNS=$1
fi

./compile.sh  # Компилируем решения

TESTS_FILE_PATH="tests.txt"
PROGRAMS_FILE_PATH="programs.txt"
INDICATOR_TYPE="2"  # 0 - Отключен, 1 - Время, 2 - "Прогрессбар"

mapfile -t tests < "$TESTS_FILE_PATH"
mapfile -t programs < "$PROGRAMS_FILE_PATH"

echo "Будем проводить по $NUM_OF_RUNS тестов"

for program in "${programs[@]}"; do

    if [ ! -e "$program" ]; then
        echo "Файл $program не найден"
        continue
    fi

    echo "Тестирование программы: $program"
    echo "=================================================="

    # Цикл по каждому входному параметру
    for test_value in "${tests[@]}"; do
        echo "Тестирование при N=$test_value"
        total_time=0
        max_time=-1
        min_time=-1
        success=1

        if [ "$INDICATOR_TYPE" -eq 1 ]; then
            echo -n "Solve times "
        elif [ "$INDICATOR_TYPE" -eq 2 ]; then
            echo -n "progress: "
        fi

        for ((i=1; i<=$NUM_OF_RUNS; i++)); do           
            { time { echo $test_value | ./$program > out.tmp; }; } 2> time.tmp

            if [ $? -ne 0 ]; then
                echo "При выполнении программы произошла ошибка"
                sed '$d' time.tmp | sed '$d' | sed '$d'
                success=0
                break
            fi

            solve_time=$(head -n 1 "out.tmp")

            if [ "$INDICATOR_TYPE" -eq 1 ]; then
                echo -n  $solve_time " "
            elif [ "$INDICATOR_TYPE" -eq 2 ]; then
                echo -n "|"
            fi

            if [ "$solve_time" -gt "$max_time" ]; then
                max_time=$solve_time
            fi
            
            if [ "$min_time" -lt 0 ] || [ "$solve_time" -lt "$min_time" ]; then
                min_time=$solve_time
            fi
            
 
            total_time=$(( total_time + solve_time ))
        done

        if [ "$success" -eq 0 ]; then
          continue
        fi

        avg_time=$(expr $total_time / $NUM_OF_RUNS)
        echo -e "\n"
        echo "MIN solve time: " $( ./time2human.sh $min_time ) " ($min_time mks)" 
        echo "MAX solve time: " $( ./time2human.sh $max_time ) " ($max_time mks)"
        echo "AVG solve time: " $( ./time2human.sh $avg_time ) " ($avg_time mks)"
        echo "error: " $(tail -n 1 "out.tmp")

        if [ $(wc -l < "out.tmp") -gt 2 ]; then
            sed '1d;$d' "out.tmp"
        fi

        cat time.tmp | grep 'real'
        echo -e "\n\n"
    done
done

rm *.tmp
