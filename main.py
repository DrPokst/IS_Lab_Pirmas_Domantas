"""
Sukurkite paprasto klasifikatoriaus (Perceptrono) išėjimui apskaičiuoti skirtą programą.
Klasifikatorius turi skirstyti objektus į dvi klases, pagal du požymius.
Išėjimo skaičiavimas atliekamas pagal formulę: y = 1, kai x1*w1 + x2*w2 + b > 0; y = -1, kai x1*w1 + x2*w2 + b <= 0;
čia w1, w2 ir b parametrai, kurie turi būti sugeneruojami naudojant atsitiktinių skaičių generatorių
(MATLAB pvz.: w1 = randn(1);) pirmąją programos veikimo iteraciją ir vėliau atnaujinami mokymo algoritmu;
x1 ir x2 yra objektų požymiai, apskaičiuoti Matlab funkcijomis, esančiomis paruoštame kodo ruošinyje arba
Data.txt faile (kiekvienoje eilutėje yra toks duomenų formatas: požymis1, požymis2, norimas_atsakymas),
jei ketinate naudoti ne Matlab..

"""
import numpy as np
from io import StringIO

# read Data file
text_file = open("Data.txt", "r")
c = StringIO(u"0 1\n2 3")
x1, x2, T = np.loadtxt(text_file, delimiter=',', usecols=(0, 1, 2), unpack=True)
x11 = np.zeros(6, )
x22 = np.zeros(6, )

# Sugeneruojami random koeficientai w1 w2 b

w1 = np.random.rand(1, )
w2 = np.random.rand(1, )
b = np.random.rand(1, )  # arba galima w = np.zeros((1,))

x11 = np.array([0.14115, 0.31565, 0.46111, 0.18474, 0.098166])
x22 = np.array([0.83535, 0.83101, 0.82518, 0.6279, 0.79092])
TT = np.array([1, 1, 1, -1, -1])



def OutputandLossCalculate(w1, w2, b, x1, x2, t):
    """
    FUunkcija paskaičiuoja parceptrono išėjimą

    Arguments:
        w1 -- x1 įėjimo svoris
        w2 -- x2 įėjimo svoris
        b -- parceptrono svoris
        x1 -- x1 iejimo duomenys
        x2 -- x2 iejimo duomenys
        T -- tikrasis atsakymas

    Returns:
        v -- išejimo reiksmė
        y -- priskirta isejimo reiksme
        e -- paskaičiuota klaida
   """
    length = len(x1)

    v = x1 * w1 + x2 * w2 + b

    y = np.zeros(length, )

    for x in range(length):
        if v[x] > 0:
            y[x] = 1
        else:
            y[x] = -1

    # Randame e
    e = t - y

    return v, y, e


def optimize(w1, w2, b, x1, x2, e, learning_rate):
    """
    This function optimizes w and b by running a gradient descent algorithm
   """
    d = 0
    length = len(x1)
    for i in range(length):
        if e[i] != 0:
            w1 = w1 + learning_rate * e[i] * x1[i]
            w2 = w2 + learning_rate * e[i] * x2[i]
            b = b + learning_rate * e[i]

            d += 1

    print("W1 W2 b skaiciavo :", d)

    return w1, w2, b

print("pradzioje W1",w1)

v, y, e = OutputandLossCalculate(w1, w2, b, x11, x22, TT)


learning_rate = 0.025
iteriation_nr = 0


while np.sum(np.abs(e)) != 0:

    w1, w2, b = optimize(w1, w2, b, x11, x22, e, learning_rate)
    v, y, e = OutputandLossCalculate(w1, w2, b, x11, x22, TT)
    iteriation_nr += 1


print("Reikalingu iteraciju skaičius ", iteriation_nr)

print("Svoris w1 ", w1)
print("Svoris w2 ", w2)
print("b: ", b)

v_g, y_g, e_g = OutputandLossCalculate(w1, w2, b, x1, x2, T)

print("Error atvaizdavimas ", e_g)
