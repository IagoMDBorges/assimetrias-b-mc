# monte carlo: repete o experimento M vezes, guarda as contagens e tira as métricas
import math
import os
import numpy as np
from src.geracao import medias_esperadas, gerar_contagens, A_total_esperado
from src.estimadores import soma, media, desvio_padrao, assimetria, sigma_A, pull


def garantir_pasta(nome):
    if not os.path.exists(nome):
        os.makedirs(nome)

def nome_arquivo(cenario, N):
    return "resultados/contagens_" + cenario + "_N" + str(N) + ".npz"

def rodar_configuracao(N, pesos, A, M, rng):
    mu_m, mu_p = medias_esperadas(N, pesos, A)
    todos_m = []
    todos_p = []
    for r in range(M):
        n_m, n_p = gerar_contagens(mu_m, mu_p, rng)
        todos_m.append(n_m)
        todos_p.append(n_p)
    return todos_m, todos_p

def salvar_contagens(caminho, todos_m, todos_p):
    garantir_pasta("resultados")
    np.savez_compressed(caminho, n_menos=np.array(todos_m), n_mais=np.array(todos_p))

def carregar_contagens(caminho):
    d = np.load(caminho)
    return d["n_menos"].tolist(), d["n_mais"].tolist()

def metricas_assimetria(todos_m, todos_p, A_reg, pesos):
    # uma entrada por região e a última é o total
    M = len(todos_m)
    k = len(pesos)
    A_tot = A_total_esperado(pesos, A_reg)
    saida = []
    for i in range(k + 1):
        if i < k:
            A_verd = A_reg[i]
        else:
            A_verd = A_tot
        ests = []
        pulls = []
        Ts = []
        for r in range(M):
            if i < k:
                nm = todos_m[r][i]
                np_ = todos_p[r][i]
            else:
                nm = soma(todos_m[r])
                np_ = soma(todos_p[r])
            a = assimetria(nm, np_)
            if a is None:
                continue
            T = nm + np_
            ests.append(a)
            Ts.append(T)
            p_val = pull(a, A_verd, T)
            if p_val is not None:
                pulls.append(p_val)
        difs = []
        for a in ests:
            difs.append(a - A_verd)
        desvio = desvio_padrao(ests)
        previsto = sigma_A(A_verd, media(Ts))
        saida.append({
            "A_verdadeiro": A_verd,
            "vies": media(difs),
            "erro_vies": desvio / math.sqrt(len(ests)),
            "desvio": desvio,
            "previsto": previsto,
            "razao": desvio / previsto,
            "pull_media": media(pulls),
            "pull_desvio": desvio_padrao(pulls),
        })
    return saida
