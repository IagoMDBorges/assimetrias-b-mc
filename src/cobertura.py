# cobertura dos intervalos: a soma exata e a contagem do monte carlo
import math
from scipy import stats
from src.estimadores import soma
from src.intervalos import intervalo_p
from src.geracao import A_total_esperado


def cobertura_exata(n, p, metodo):
    # soma, sobre k, da pmf vezes (o intervalo de k contém p?)
    centro = n * p
    desvio = math.sqrt(n * p * (1 - p))
    k_ini = int(centro - 7 * desvio)
    k_fim = int(centro + 7 * desvio) + 1
    if k_ini < 0:
        k_ini = 0
    if k_fim > n:
        k_fim = n
    total = 0.0
    for k in range(k_ini, k_fim + 1):
        inf, sup = intervalo_p(k, n, metodo)
        if inf <= p and p <= sup:
            total = total + stats.binom.pmf(k, n, p)
    return float(total)


def cobertura_media(lam, p, metodo, max_pontos=21):
    # T também flutua (poisson), então tiro a média da cobertura sobre T
    desvio = math.sqrt(lam)
    t_ini = int(lam - 5 * desvio)
    if t_ini < 1:
        t_ini = 1
    t_fim = int(lam + 5 * desvio) + 1
    if lam <= 200:
        passo = 1
    else:
        passo = (t_fim - t_ini) // max_pontos + 1
    soma_cob = 0.0
    soma_peso = 0.0
    t = t_ini
    while t <= t_fim:
        peso = stats.poisson.pmf(t, lam)
        soma_cob = soma_cob + peso * cobertura_exata(t, p, metodo)
        soma_peso = soma_peso + peso
        t = t + passo
    return float(soma_cob / soma_peso)


def cobertura_mc(todos_m, todos_p, A_reg, pesos, metodo):
    # devolve (cobertura das regiões juntas, nº de tentativas, cobertura do total)
    M = len(todos_m)
    k = len(pesos)
    A_tot = A_total_esperado(pesos, A_reg)
    ac_reg = 0
    tent_reg = 0
    ac_tot = 0
    for r in range(M):
        for i in range(k):
            T = todos_m[r][i] + todos_p[r][i]
            if T == 0:
                continue
            p_verd = (1 + A_reg[i]) / 2
            inf, sup = intervalo_p(todos_m[r][i], T, metodo)
            tent_reg = tent_reg + 1
            if inf <= p_verd and p_verd <= sup:
                ac_reg = ac_reg + 1
        nm = soma(todos_m[r])
        np_ = soma(todos_p[r])
        p_verd = (1 + A_tot) / 2
        inf, sup = intervalo_p(nm, nm + np_, metodo)
        if inf <= p_verd and p_verd <= sup:
            ac_tot = ac_tot + 1
    return ac_reg / tent_reg, tent_reg, ac_tot / M


def cobertura_exata_regioes(N, pesos, A_reg, metodo, guardado):
    # média da cobertura exata nas regiões; "guardado" evita refazer conta igual
    total = 0
    for i in range(len(pesos)):
        lam = N * pesos[i]
        p = (1 + A_reg[i]) / 2
        chave = (metodo, round(lam, 3), round(p, 6))
        if chave not in guardado:
            guardado[chave] = cobertura_media(lam, p, metodo)
        total = total + guardado[chave]
    return total / len(pesos)
