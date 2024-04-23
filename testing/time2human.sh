microseconds=$1

# Разделяем микросекунды на секунды и микросекунды
seconds=$((microseconds / 1000000))
micro=$((microseconds % 1000000))

# Разделяем секунды на минуты и секунды
minutes=$((seconds / 60))
seconds=$((seconds % 60))

# Отделяем миллисекунды от микросекунд
milliseconds=$((micro / 1000))

# Добавляем нули спереди, если необходимо
printf -v seconds "%02d" "$seconds"
printf -v milliseconds "%03d" "$milliseconds"

# Выводим результат
echo "${minutes}m${seconds}.${milliseconds}s"


