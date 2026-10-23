# sorteio das contagens (poisson) e leitura do config
import json
import numpy as np


def ler_config(caminho="config.json"):
    arq = open(caminho)
    cfg = json.load(arq)
    arq.close()
    # confere o básico antes de sortear qualquer coisa
    total = 0
    for w in cfg["pesos"]:
        total = total + w
    if abs(total - 1) > 1e-9:
        raise ValueError("os pesos não somam 1")
    for nome in cfg["cenarios"]:
        if len(cfg["cenarios"][nome]) != cfg["k"]:
            raise ValueError("tamanho errado no cenário " + nome)
        for A in cfg["cenarios"][nome]:
            if abs(A) > 1:
                raise ValueError("assimetria maior que 1 no cenário " + nome)
    return cfg

def criar_geradores(semente_mae, quantos):
    # cada configuração usa o seu gerador, todos nascem da mesma semente-mãe
    filhos = np.random.SeedSequence(semente_mae).spawn(quantos)
    lista = []
    for f in filhos:
        lista.append(np.random.default_rng(f))
    return lista

def medias_esperadas(N, pesos, A):
    mu_m = []
    mu_p = []
    for i in range(len(pesos)):
        metade = N * pesos[i] / 2
        mu_m.append(metade * (1 + A[i]))
        mu_p.append(metade * (1 - A[i]))
    return mu_m, mu_p


def A_total_esperado(pesos, A):
    t = 0
    for i in range(len(pesos)):
        t = t + pesos[i] * A[i]
    return round(t, 12)

def gerar_contagens(mu_m, mu_p, rng):
    n_m = []
    n_p = []
    for i in range(len(mu_m)):
        n_m.append(int(rng.poisson(mu_m[i])))
        n_p.append(int(rng.poisson(mu_p[i])))
    return n_m, n_p
