import sympy as sym


EI, ll, F, a = sym.symbols('EI, l, F, a')
l = a * ll
# Массив коэффициентов у EI, l для первого и второго конечного элемента
args = [2, 2, 1, 1]
# Массив для ограничений в узлах. 0 - без приколов, 1 - заделка, 2 - шарнир односвязный, 3 - шарнир двусвязный,
# 4 - паз, 5 - поршень
ligma = [5, 2, 1]
# Массив сил на узлах.
forces = sym.Matrix([-F, 0, 0, -F * ll, 0, 0])


# Получает граничные условия из ограничений в узлах
def Deligmitisator(lig):
    deligma = []
    for i in lig:
        if i == 0:
            deligma.append(None)
            deligma.append(None)
        if i == 1:
            deligma.append(0)
            deligma.append(0)
        if i == 2:
            deligma.append(0)
            deligma.append(None)
        if i == 3:
            deligma.append(0)
            deligma.append(None)
        if i == 4:
            deligma.append(0)
            deligma.append(0)
        if i == 5:
            deligma.append(None)
            deligma.append(0)
    return deligma


# Подсчет матрицы жесткости для конечного элемента
def Matrix(EI1, l1):
    k = sym.Matrix([[12 * EI1 * EI / (l1 * l) ** 3, -6 * EI1 * EI / (l1 * l) ** 2, -12 * EI1 * EI / (l1 * l) ** 3,
                     -6 * EI1 * EI / (l1 * l) ** 2],
                    [-6 * EI1 * EI / (l1 * l) ** 2, 4 * EI1 * EI / (l1 * l), 6 * EI1 * EI / (l1 * l) ** 2,
                     2 * EI1 * EI / (l1 * l)],
                    [-12 * EI1 * EI / (l1 * l) ** 3, 6 * EI1 * EI / (l1 * l) ** 2, 12 * EI1 * EI / (l1 * l) ** 3,
                     6 * EI1 * EI / (l1 * l) ** 2],
                    [-6 * EI1 * EI / (l1 * l) ** 2, 2 * EI1 * EI / (l1 * l), 6 * EI1 * EI / (l1 * l) ** 2,
                     4 * EI1 * EI / (l1 * l)]])
    return k / (EI / l ** 3)


# Подсчет ансамбля
def Revansamble(k1, k2):
    ans = sym.zeros(6, 6)
    for i in range(4):
        for j in range(4):
            ans[i, j] = k1[i, j]
    for i in range(2, 6):
        for j in range(2, 6):
            ans[i, j] += k2[i - 2, j - 2]
    return ans


# Шикарный вывод матрицы
def Mprint(k):
    for i in range(len(k.row(0))):
        for j in range(k.shape[0]):
            text = k[i, j]
            print(str(text) + (15 - len(str(text))) * ' ' + '|', end=' ')
        print()
    print()


# Дезинтегратор ансамблей
def Desintegrator(k, chokma, force):
    for i in range(len(chokma)):
        if chokma[i] == 0:
            force[i] = 0
            for j in range(len(k.row(i))):
                if j != i:
                    k[i, j] = 0
                    k[j, i] = 0
                else:
                    k[i, j] = 1


print("Перед каждой матрицей впишите коэффициент EI / l^3\n")
k_1 = Matrix(*args[:2])
print("Матрица k1")
Mprint(k_1)
k_2 = Matrix(*args[2:])
print("Матрица k2")
Mprint(k_2)

print("Ансамбль без учета граничных условий")
ansamble = Revansamble(k_1, k_2)
Mprint(ansamble)

suckma = Deligmitisator(ligma)
Desintegrator(ansamble, suckma, forces)
congratulations = ['w1', 'o1', 'w2', 'o2', 'w3', 'o3']
well_done = ['f11', 'f12', 'f21', 'f22', 'f31', 'f32']

print("Граничные условия")
for i in range(len(suckma)):
    if suckma[i] == 0:
        print(f"{congratulations[i]} = {suckma[i]}")
    if forces[i] != 0:
        print(f"{well_done[i]} = {forces[i]}")

print("\nУльтра Ансамбль")
Mprint(ansamble)
u = (ansamble ** -1 * forces / (EI / l ** 3)).col(0)
print("Вектор u:")
for i in range(len(u)):
    print(f"{congratulations[i]} = {u[i]}")
