# estimadores da assimetria, incertezas e umas funções de apoio
import math


def soma(lista):
    t = 0
    for x in lista:
        t = t + x
    return t

def media(lista):
    return soma(lista) / len(lista)

def desvio_padrao(lista):
    m = media(lista)
    acc = 0
    for x in lista:
        acc = acc + (x - m) ** 2
    return math.sqrt(acc / (len(lista) - 1))

def inclinacao(xs, ys):
    mx = media(xs)
    my = media(ys)
    num = 0
    den = 0
    for i in range(len(xs)):
        num = num + (xs[i] - mx) * (ys[i] - my)
        den = den + (xs[i] - mx) ** 2
    return num / den


def assimetria(n_menos, n_mais):
    T = n_menos + n_mais
    if T == 0:
        return None     
    return (n_menos - n_mais) / T

def assimetria_total(lista_menos, lista_mais):
    return assimetria(soma(lista_menos), soma(lista_mais))


def assimetria_ponderada(lista_menos, lista_mais):
    # o mesmo total de antes, agora como média das locais ponderada por T_i
    num = 0
    den = 0
    for i in range(len(lista_menos)):
        T = lista_menos[i] + lista_mais[i]
        a = assimetria(lista_menos[i], lista_mais[i])
        if a is not None:
            num = num + T * a
            den = den + T
    if den == 0:
        return None
    return num / den


def media_simples(lista_menos, lista_mais):
    # média sem peso nenhum, só pra comparar. não serve como estimador do total
    vals = []
    for i in range(len(lista_menos)):
        a = assimetria(lista_menos[i], lista_mais[i])
        if a is not None:
            vals.append(a)
    return media(vals)

def sigma_A(A, T):
    return math.sqrt((1 - A * A) / T)

def pull(A_est, A_verd, T):
    if A_est is None or abs(A_est) >= 1:
        return None      # sigma estimado dá zero
    return (A_est - A_verd) / sigma_A(A_est, T)

# ---- máxima verossimilhança da binomial ----
def log_verossimilhanca(p, k, n):
    return k * math.log(p) + (n - k) * math.log(1 - p)

def emv_binomial(k, n):
    # busca ternária do máximo; tem que dar k/n
    baixo = 1e-9
    alto = 1 - 1e-9
    for passo in range(200):
        a = baixo + (alto - baixo) / 3
        b = alto - (alto - baixo) / 3
        if log_verossimilhanca(a, k, n) < log_verossimilhanca(b, k, n):
            baixo = a
        else:
            alto = b
    return (baixo + alto) / 2

def informacao_fisher(p, n):
    return n / (p * (1 - p))

def variancia_numerica(k, n, h=1e-4):
    # 1 / (-d2 lnL / dp2) em p = k/n, com diferença finita
    p = k / n
    d2 = (log_verossimilhanca(p + h, k, n) - 2 * log_verossimilhanca(p, k, n) + log_verossimilhanca(p - h, k, n)) / (h * h)
    return -1 / d2
