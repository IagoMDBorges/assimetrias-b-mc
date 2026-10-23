# testes de hipótese: valor-p exato, falso positivo, poder, várias regiões e qui-quadrado
import math
from scipy import stats
from src.estimadores import soma
from src.geracao import medias_esperadas, gerar_contagens


def pvalor_exato(k, n):
    # H0: A = 0, ou seja p = 1/2. binomial bilateral
    menor = min(k, n - k)
    p = 2 * stats.binom.cdf(menor, n, 0.5)
    if p > 1:
        p = 1.0
    return float(p)


def taxa_rejeicao_exata(n, A, alfa):
    # soma da pmf dos k que levam a rejeitar H0: com A=0 é o falso positivo, com A!=0 é o poder
    p = (1 + A) / 2
    centro = n * p
    desvio = math.sqrt(n * p * (1 - p))
    k_ini = int(centro - 8 * desvio)
    k_fim = int(centro + 8 * desvio) + 1
    if k_ini < 0:
        k_ini = 0
    if k_fim > n:
        k_fim = n
    total = 0.0
    for k in range(k_ini, k_fim + 1):
        if pvalor_exato(k, n) <= alfa + 1e-12:
            total = total + stats.binom.pmf(k, n, p)
    return float(total)


def fp_exato(n, alfa):
    return taxa_rejeicao_exata(n, 0.0, alfa)


def poder_exato(n, A, alfa):
    return taxa_rejeicao_exata(n, A, alfa)


def poder_normal(A, T):
    desvio = abs(A) * math.sqrt(T)
    lado_direito = stats.norm.cdf(desvio - 1.96)
    lado_esquerdo = stats.norm.cdf(-desvio - 1.96)
    return float(lado_direito + lado_esquerdo)


def z_regiao(n_m, n_p):
    T = n_m + n_p
    if T == 0:
        return None
    return (n_m - n_p) / math.sqrt(T)


def chi2_zero(lista_m, lista_p):
    # soma dos z^2: testa A_i = 0 em todas as regiões
    total = 0.0
    usadas = 0
    for i in range(len(lista_m)):
        z = z_regiao(lista_m[i], lista_p[i])
        if z is not None:
            total = total + z * z
            usadas = usadas + 1
    return total, usadas


def chi2_hom(lista_m, lista_p):
    # testa se a assimetria é a mesma em todas as regiões (S_i com o alfa global)
    alfa = soma(lista_m) / soma(lista_p)
    total = 0.0
    usadas = 0
    for i in range(len(lista_m)):
        T = lista_m[i] + lista_p[i]
        if T == 0:
            continue
        S = (lista_m[i] - alfa * lista_p[i]) / math.sqrt(alfa * T)
        total = total + S * S
        usadas = usadas + 1
    return total, usadas - 1


def max_abs_z(lista_m, lista_p):
    maior = 0.0
    for i in range(len(lista_m)):
        z = z_regiao(lista_m[i], lista_p[i])
        if z is not None and abs(z) > maior:
            maior = abs(z)
    return maior


def pvalores_regioes(lista_m, lista_p):
    ps = []
    for i in range(len(lista_m)):
        T = lista_m[i] + lista_p[i]
        if T > 0:
            ps.append(pvalor_exato(lista_m[i], T))
    return ps


def bonferroni_rejeita(ps, alfa):
    return min(ps) <= alfa / len(ps)


def p_global(valor_obs, nulos):
    # p global calibrado por monte carlo: (1 + #nulos >= obs) / (1 + M0)
    cont = 0
    for v in nulos:
        if v >= valor_obs:
            cont = cont + 1
    return (1 + cont) / (1 + len(nulos))


def gerar_nulos(N, pesos, M0, rng):
    # experimentos com A = 0 em tudo, pra calibrar as estatísticas globais
    mu_m, mu_p = medias_esperadas(N, pesos, [0.0] * len(pesos))
    maxz = []
    chi0 = []
    chih = []
    for r in range(M0):
        n_m, n_p = gerar_contagens(mu_m, mu_p, rng)
        maxz.append(max_abs_z(n_m, n_p))
        chi0.append(chi2_zero(n_m, n_p)[0])
        chih.append(chi2_hom(n_m, n_p)[0])
    return maxz, chi0, chih
