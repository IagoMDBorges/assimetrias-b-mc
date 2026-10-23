# eficiência diferente entre B- e B+ imita uma assimetria
from src.estimadores import soma


def A_obs(A, a):
    # a = (eps- - eps+) / (eps- + eps+)
    return (A + a) / (1 + A * a)


def corrigir(A_observado, a):
    return (A_observado - a) / (1 - A_observado * a)


def medias_com_eficiencia(N, pesos, A, a):
    mu_m = []
    mu_p = []
    for i in range(len(pesos)):
        metade = N * pesos[i] / 2
        mu_m.append(metade * (1 + A[i]) * (1 + a))
        mu_p.append(metade * (1 - A[i]) * (1 - a))
    return mu_m, mu_p


def A_obs_total_esperado(N, pesos, A, a):
    mu_m, mu_p = medias_com_eficiencia(N, pesos, A, a)
    return (soma(mu_m) - soma(mu_p)) / (soma(mu_m) + soma(mu_p))
