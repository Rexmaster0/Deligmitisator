import os
import subprocess
import sympy as sym
import shutil

EI, l, F, a = sym.symbols('EI, l, F, a')
types = {"ничего": 0, "заделка": 1, "шарнир односвязный": 2, "шарнир двусвязный": 3, "паз": 4, "поршень": 5}

# Массив коэффициентов у EI, l для первого и второго конечного элемента
args = [2, 2 * a, 1, a]
# Массив для типов ограничений. Введите названия в строку
limits = "поршень, шарнир односвязный, заделка"
# Массив сил на узлах. Неизвестные силы заносите как нулевые
forces = [-F, 0, 0, -F*l, 0, 0]

if str(a) in ''.join(str(args)):
    flag = True
    div = (EI / (l * a) ** 3)
else:
    div = (EI / l ** 3)


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

    return k / div


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


print(f"Перед каждой матрицей впишите коэффициент {div}\n")
k_1 = Matrix(*args[:2])
print("Матрица k1")
Mprint(k_1)
k_2 = Matrix(*args[2:])
print("Матрица k2")
Mprint(k_2)

print("Ансамбль без учета граничных условий")
ansamble = Revansamble(k_1, k_2)
Mprint(ansamble)
print("Граничные условия")
ligma = [types[i] for i in limits.split(', ')]
forces = sym.Matrix(forces)
suckma = Deligmitisator(ligma)
congratulations = [r'$\omega_{1}$', r'$\theta_{1}$', r'$\omega_{2}$', r'$\theta_{2}$', r'$\omega_{3}$', r'$\theta_{3}$']
well_done = ['$f_{11}$', '$f_{12}$', '$f_{21}$', '$f_{22}$', '$f_{31}$', '$f_{32}$']
restrict = ''
for i in range(len(suckma)):
    if suckma[i] == 0:
        print(f"{congratulations[i]} = {suckma[i]}    {well_done[i]} = ?")
        restrict += fr"{congratulations[i]} = {suckma[i]} \quad {well_done[i]} = ?\\"
    else:
        print(f"{congratulations[i]} = ?    {well_done[i]} = {sym.latex(forces[i])}")
        restrict += fr"{congratulations[i]} = ? \quad {well_done[i]} = {sym.latex(forces[i])} \\"
ansamblecopy = ansamble.copy()
print("\nУльтра Ансамбль")
Desintegrator(ansamble, suckma, forces)
Mprint(ansamble)
u = (ansamble ** -1 * forces / div).col(0)

print("Вектор f:")
for i in range(len(forces)):
    print(f"{well_done[i]} = {forces[i]}")

print("\nВектор u:")
for i in range(len(u)):
    print(f"{congratulations[i]} = {u[i]}")

w1, o1, w2, o2, w3, o3 = sym.symbols(r'$\omega_{1}$, $\theta_{1}$, $\omega_{2}$, $\theta_{2}$, $\omega_{3}$, '
                                     r'$\theta_{3}$')
linguine = sym.Matrix([w1, o1, w2, o2, w3, o3])
for i in range(len(suckma)):
    if suckma[i] == 0:
        linguine[i] = 0

with open('ligma.tex', 'r', encoding='utf-8') as file:
    data = file.readlines()
data[23] = sym.latex(div)
data[24] = sym.latex(k_1)
data[27] = sym.latex(div)
data[28] = sym.latex(k_2)
data[32] = sym.latex(div) + sym.latex(ansamblecopy)
data[35] = restrict
data[39] = sym.latex(div) + sym.latex(ansamblecopy)
data[59] = sym.latex(div) + sym.latex(ansamble)
data[60] = sym.latex(linguine, mat_delim='(').replace('$', '')
data[62] = sym.latex(forces, mat_delim='(')
count = 0
for i in range(len(forces)):
    if linguine[i] != 0:
        data.insert(68 + count, (sym.latex((div * ansamble.row(i) * linguine), mat_delim='') +
                                 fr" = {sym.latex(forces[i])} \\ ").replace('$', ''))
        count += 1
data[68 + count] = data[68 + count][:-3]
for i in range(len(u)):
    if u[i] != 0:
        data.insert(70 + count, fr"{sym.latex(linguine[i])} = {sym.latex(u[i])} \quad \quad".replace('$', ''))
        count += 1

data[81 + count] = " = " + sym.latex(u, mat_delim='(')
file.close()
path = os.getcwd()
with open('rk2.tex', 'w', encoding='utf-8') as file:
    for i in data:
        if i[-1] == "\n":
            file.write(i)
        else:
            file.write(i + "\n")

print("\nЕсли у вас не Windows, то для переноса rk2.tex в pdf вам надо закинуть его на этот сайт:"
      "\nhttps://products.groupdocs.app/conversion/tex-to-pdf")
os.chdir(fr"{path}\texlive_minimal")
subprocess.run(['pdflatex.exe', '-interaction=nonstopmode', fr'{path}\rk2.tex'], stdout=subprocess.DEVNULL,
               stderr=subprocess.DEVNULL)
shutil.move("rk2.pdf", rf"{path}\rk2.pdf")
