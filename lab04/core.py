import random
import time


class LinearCongruentialGenerator:
    def __init__(self, seed, a=16807, c=12345, m=(2 ** 31 - 1)):
        self.state = seed
        self.a = a
        self.c = c
        self.m = m

    def random(self):
        self.state = (self.a * self.state + self.c) % self.m
        return self.state / self.m


def calculate_mean(data):
    return sum(data) / len(data)

def calculate_variance(data, mean):
    n = len(data)
    return sum((x - mean) ** 2 for x in data) / (n - 1)


def main():
    print("╔" + "═"*78 + "╗")
    print("║" + "ЛАБОРАТОРНАЯ РАБОТА: БАЗОВЫЙ ДАТЧИК СЛУЧАЙНЫХ ЧИСЕЛ".center(78) + "║")
    print("╚" + "═"*78 + "╝\n")

    print("► Настройка параметров")
    try:
        seed_input = input("  Введите начальное значение (seed) [Enter для генерации по времени]: ")
        seed = int(seed_input) if seed_input.strip() else int(time.time() * 1000)

        sample_size_input = input("  Введите размер выборки (sample_size) [Enter для генерации по времени]: ")
        sample_size = int(sample_size_input) if sample_size_input.strip() else 100_000

    except ValueError:
        print("  [!] Некорректный ввод. Использовано значение по времени.")
        seed = int(time.time() * 1000)
        sample_size = 100_000

    print(f"  Размер выборки: {sample_size:,} значений\n".replace(',', ' '))

    custom_rng = LinearCongruentialGenerator(seed, c=0)
    random.seed(seed)

    print("► Генерация чисел...")
    custom_sample = [custom_rng.random() for _ in range(sample_size)]
    builtin_sample = [random.random() for _ in range(sample_size)]

    print("► Вычисление статистических характеристик...\n")
    custom_mean = calculate_mean(custom_sample)
    custom_var = calculate_variance(custom_sample, custom_mean)

    builtin_mean = calculate_mean(builtin_sample)
    builtin_var = calculate_variance(builtin_sample, builtin_mean)

    theory_mean = 0.5
    theory_var = 1 / 12

    print("╔" + "═"*78 + "╗")
    print("║" + "РЕЗУЛЬТАТЫ ЭКСПЕРИМЕНТА".center(78) + "║")
    print("╠" + "═"*26 + "╦" + "═"*25 + "╦" + "═"*25 + "╣")
    print("║" + " Генератор".ljust(26) + "║" + " Математическое ожидание ".center(25) + "║" + " Дисперсия ".center(25) + "║")
    print("╠" + "═"*26 + "╬" + "═"*25 + "╬" + "═"*25 + "╣")
    print(f"║ Теоретические значения   ║ {theory_mean:^23.6f} ║ {theory_var:^23.6f} ║")
    print("╠" + "═"*26 + "╬" + "═"*25 + "╬" + "═"*25 + "╣")
    print(f"║ Базовый (LCG)            ║ {custom_mean:^23.6f} ║ {custom_var:^23.6f} ║")
    print(f"║ Погрешность (отн. теор.) ║ {abs(custom_mean - theory_mean):^23.6f} ║ {abs(custom_var - theory_var):^23.6f} ║")
    print("╠" + "═"*26 + "╬" + "═"*25 + "╬" + "═"*25 + "╣")
    print(f"║ Встроенный (Python)      ║ {builtin_mean:^23.6f} ║ {builtin_var:^23.6f} ║")
    print(f"║ Погрешность (отн. теор.) ║ {abs(builtin_mean - theory_mean):^23.6f} ║ {abs(builtin_var - theory_var):^23.6f} ║")
    print("╚" + "═"*26 + "╩" + "═"*25 + "╩" + "═"*25 + "╝\n")

if __name__ == "__main__":
    main()
