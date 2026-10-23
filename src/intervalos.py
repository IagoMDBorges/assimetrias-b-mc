# intervalos de confiança para p = N-/T. no fim viram intervalos de A = 2p - 1
import math
from scipy import stats
from src.estimadores import log_verossimilhanca

Z95 = 1.959963984540054


def wald(k, n, z=Z95):
    p = k / n
    meia = z * math.sqrt(p * (1 - p) / n)
    return p - meia, p + meia


def wilson(k, n, z=Z95):
    p = k / n
    z2 = z * z
    centro = (p + z2 / (2 * n)) / (1 + z2 / n)
    meia = z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / (1 + z2 / n)
    return centro - meia, centro + meia


def clopper_pearson(k, n, nivel=0.95):
    a = 1 - nivel
    if k == 0:
        inf = 0.0
    else:
        inf = stats.beta.ppf(a / 2, k, n - k + 1)
    if k == n:
        sup = 1.0
    else:
        sup = stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return float(inf), float(sup)


def razao_verossimilhanca(k, n, nivel=0.95):
    # região onde 2*(lnL(p_hat) - lnL(p)) <= quantil do qui-quadrado com 1 g.l.
    corte = stats.chi2.ppf(nivel, 1) / 2
    if k == 0:
        return 0.0, 1 - math.exp(-corte / n)
    if k == n:
        return math.exp(-corte / n), 1.0
    p_hat = k / n
    topo = log_verossimilhanca(p_hat, k, n)
    a = 1e-12
    b = p_hat
    for passo in range(80):
        meio = (a + b) / 2
        if log_verossimilhanca(meio, k, n) < topo - corte:
            a = meio
        else:
            b = meio
    inf = (a + b) / 2
    a = p_hat
    b = 1 - 1e-12
    for passo in range(80):
        meio = (a + b) / 2
        if log_verossimilhanca(meio, k, n) < topo - corte:
            b = meio
        else:
            a = meio
    sup = (a + b) / 2
    return inf, sup


def intervalo_p(k, n, metodo):
    if metodo == "wald":
        return wald(k, n)
    elif metodo == "wilson":
        return wilson(k, n)
    elif metodo == "clopper_pearson":
        return clopper_pearson(k, n)
    elif metodo == "verossimilhanca":
        return razao_verossimilhanca(k, n)
    else:
        raise ValueError("método desconhecido: " + metodo)


def p_para_A(p_inf, p_sup):
    return 2 * p_inf - 1, 2 * p_sup - 1


def intervalo_A(n_menos, n_mais, metodo):
    T = n_menos + n_mais
    if T == 0:
        return None
    p_inf, p_sup = intervalo_p(n_menos, T, metodo)
    A_inf, A_sup = p_para_A(p_inf, p_sup)
    if A_inf < -1:
        A_inf = -1.0
    if A_sup > 1:
        A_sup = 1.0
    return A_inf, A_sup
