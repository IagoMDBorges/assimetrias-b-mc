# roda os blocos S1 a S11 em sequência e grava tudo em resultados/ e figuras/
# leva uns 2-3 min. a parte lenta é a cobertura exata 
import csv
import math
import time
from scipy import stats

from src.geracao import ler_config, criar_geradores, medias_esperadas, A_total_esperado, gerar_contagens
from src.estimadores import (soma, media, desvio_padrao, assimetria, assimetria_total, assimetria_ponderada, media_simples, sigma_A, emv_binomial, informacao_fisher, variancia_numerica, inclinacao)
from src.intervalos import wald, wilson, clopper_pearson, intervalo_p, intervalo_A
from src.mc import (rodar_configuracao, salvar_contagens, carregar_contagens, nome_arquivo, metricas_assimetria, garantir_pasta)
from src.cobertura import cobertura_exata, cobertura_media, cobertura_mc, cobertura_exata_regioes
from src.testes import (pvalor_exato, fp_exato, poder_exato, poder_normal, chi2_zero, chi2_hom, max_abs_z, pvalores_regioes, bonferroni_rejeita, p_global, gerar_nulos)
from src.eficiencia import A_obs, corrigir, medias_com_eficiencia, A_obs_total_esperado
from src.espaco_fase import (m1, m2, m3, M_B, s23_limites, pesos_por_faixa, centros_das_faixas, perfil_gaussiano, largura_otima, significancia_esperada)
from src.dalitz import (gerar_posicoes, erro_conservacao, atribuir_cargas, contar_por_bin, S_com_alfa, p_do_chi2_hom, rodar_nulos)
from src.figuras import (fig_contagens, fig_assimetrias, fig_assimetrias_ic, fig_hist_total, fig_compensacao, fig_sigma_vs_N, fig_cobertura_curva, fig_cobertura_mc, fig_poder, fig_eficiencia, fig_largura, fig_mapa_S, fig_hist_S)

cfg = ler_config()
pesos = cfg["pesos"]
Ns = cfg["N_exp"]
Ns_escala = cfg["N_escala"]
M = cfg["M"]
alfa = cfg["alfa"]
cenarios = cfg["cenarios"]
tol = cfg["tolerancias"]
METODOS = ["wald", "wilson", "clopper_pearson"]

verif = []        # tabela de verificação (etapa 1)
verif3 = []       # etapas 2 e 3
metricas = {}     # (cenário, N) -> métricas do estimador
cob = {}          # (cenário, N, método) -> cobertura MC e exata
rej = {}          # (cenário, N) -> taxas de rejeição dos testes


def linha(lista, item, descricao, valor, criterio, ok):
    if ok:
        lista.append([item, descricao, valor, criterio, "sim"])
    else:
        lista.append([item, descricao, valor, criterio, "nao"])


def s1_geracao():
    print("S1: geração")
    mu_m, mu_p = medias_esperadas(10000, pesos, [0.05] * 10)
    print("  N=10000 e A=0.05 dá mu- =", mu_m[0], "e mu+ =", mu_p[0])
    ok = True
    for nome in cenarios:
        for N in Ns:
            a, b = medias_esperadas(N, pesos, cenarios[nome])
            if abs(soma(a) + soma(b) - N) > 1e-6:
                ok = False
    for nome, esperado in [("sem_efeito", 0.0), ("global", 0.05), ("compensacao", 0.0)]:
        if abs(A_total_esperado(pesos, cenarios[nome]) - esperado) > 1e-12:
            ok = False
    linha(verif, 4, "soma dos mu = N_exp e A_total esperado", "mu-=525 e mu+=475 no exemplo", "igualdade", ok)
    # mesma semente, mesmos números
    g1 = criar_geradores(cfg["semente_mae"], 20)
    g2 = criar_geradores(cfg["semente_mae"], 20)
    a = rodar_configuracao(10000, pesos, cenarios["global"], 50, g1[4])
    b = rodar_configuracao(10000, pesos, cenarios["global"], 50, g2[4])
    linha(verif, 16, "reprodutibilidade (mesma semente)", "dois sorteios de 50 repetições", "igualdade exata", a == b)


def s2_estimadores():
    print("S2: estimadores e os três cenários")
    ok1 = abs(assimetria(525, 475) - 0.05) < 1e-12 and abs(assimetria(475, 525) + 0.05) < 1e-12
    linha(verif, 1, "convenção de sinal", "525/475 dá +0.05, trocando dá -0.05", "igualdade", ok1)
    ok2 = assimetria(3, 1) == 0.5 and assimetria(7, 7) == 0 and assimetria(5, 0) == 1 and assimetria(0, 0) is None
    linha(verif, 2, "casos feitos à mão", "(3,1) (7,7) (5,0) (0,0)", "igualdade; (0,0) indefinido", ok2)
    gs = criar_geradores(cfg["semente_mae"], 20)
    c = 0
    for nome in cenarios:
        A = cenarios[nome]
        mu_m, mu_p = medias_esperadas(10000, pesos, A)
        n_m, n_p = gerar_contagens(mu_m, mu_p, gs[c])
        c = c + 1
        A_est = []
        inf = []
        sup = []
        for i in range(len(pesos)):
            a = assimetria(n_m[i], n_p[i])
            s = sigma_A(a, n_m[i] + n_p[i])
            A_est.append(a)
            inf.append(a - 1.96 * s)
            sup.append(a + 1.96 * s)
        meta1 = assimetria(soma(n_m[0:5]), soma(n_p[0:5]))
        meta2 = assimetria(soma(n_m[5:10]), soma(n_p[5:10]))
        print("  ", nome, "| A_total esperado", round(A_total_esperado(pesos, A), 3), "estimado", round(assimetria_total(n_m, n_p), 4), "| metades", round(meta1, 3), round(meta2, 3))
        fig_contagens(n_m, n_p, "F1 - contagens (" + nome + ")", "F1_contagens_" + nome)
        fig_assimetrias(A_est, inf, sup, A, "F2 - assimetria por região (" + nome + ")", "F2_assimetrias_" + nome)


def s3_monte_carlo():
    print("S3: monte carlo")
    gs = criar_geradores(cfg["semente_mae"], 30)
    sigmas = {}
    c = 0
    falhas_vies = 0
    falhas_razao = 0
    total_m = 0
    maior_dif = 0
    for nome in cenarios:
        A = cenarios[nome]
        sigmas[nome] = []
        for N in Ns:
            todos_m, todos_p = rodar_configuracao(N, pesos, A, M, gs[c])
            c = c + 1
            salvar_contagens(nome_arquivo(nome, N), todos_m, todos_p)
            m = metricas_assimetria(todos_m, todos_p, A, pesos)
            metricas[(nome, N)] = m
            for i in range(len(m)):
                total_m = total_m + 1
                if abs(m[i]["vies"]) > tol["sigmas_vies"] * m[i]["erro_vies"]:
                    falhas_vies = falhas_vies + 1
                if abs(m[i]["razao"] - 1) > tol["razao_desvio"]:
                    falhas_razao = falhas_razao + 1
            tot = m[len(pesos)]
            sigmas[nome].append(tot["desvio"])
            print("  ", nome, N, "| total: viés", round(tot["vies"], 5), "desvio", round(tot["desvio"], 5), "previsto", round(tot["previsto"], 5))
            if N == 10000:
                lista = []
                for r in range(M):
                    lista.append(assimetria_total(todos_m[r], todos_p[r]))
                    d = abs(lista[r] - assimetria_ponderada(todos_m[r], todos_p[r]))
                    if d > maior_dif:
                        maior_dif = d
                fig_hist_total(lista, tot["A_verdadeiro"], tot["previsto"], "F3 - A total (" + nome + ", N=10000)", "F3_hist_total_" + nome)
                if nome == "compensacao":
                    meta1 = []
                    meta2 = []
                    for r in range(M):
                        meta1.append(assimetria(soma(todos_m[r][0:5]), soma(todos_p[r][0:5])))
                        meta2.append(assimetria(soma(todos_m[r][5:10]), soma(todos_p[r][5:10])))
                    fig_compensacao(lista, meta1, meta2, "F4 - compensação local (N=10000)", "F4_compensacao")

    # rodo também N=3000 pra reta em log N ter 4 pontos não simétricos
    sigmas_escala = {}
    c_extra = 23
    for nome in cenarios:
        A = cenarios[nome]
        tm_3000, tp_3000 = rodar_configuracao(3000, pesos, A, M, gs[c_extra])
        c_extra = c_extra + 1
        m_3000 = metricas_assimetria(tm_3000, tp_3000, A, pesos)[len(pesos)]
        sigmas_escala[nome] = [sigmas[nome][0], m_3000["desvio"], sigmas[nome][1], sigmas[nome][2]]
    fig_sigma_vs_N(Ns_escala, sigmas_escala, "F5_sigma_vs_N")

    linha(verif, 3, "total = média ponderada por T_i", "maior diferença " + str(maior_dif), "< 1e-12", maior_dif < 1e-12)
    # média e variância das contagens (global, N=1e4)
    todos_m, todos_p = carregar_contagens(nome_arquivo("global", 10000))
    mu_m, mu_p = medias_esperadas(10000, pesos, cenarios["global"])
    falhas = 0
    for i in range(len(pesos)):
        col_m = []
        col_p = []
        for r in range(M):
            col_m.append(todos_m[r][i])
            col_p.append(todos_p[r][i])
        for valores, mu in [(col_m, mu_m[i]), (col_p, mu_p[i])]:
            if abs(media(valores) - mu) > 3 * math.sqrt(mu / M):
                falhas = falhas + 1
            if abs(desvio_padrao(valores) ** 2 / media(valores) - 1) > 0.10:
                falhas = falhas + 1
    linha(verif, 5, "média e variância das contagens (global, N=1e4)", "fora do critério: " + str(falhas) + " de 40", "<= 2", falhas <= 2)
    linha(verif, 6, "viés de A_i e do total", "fora de 3 sigma: " + str(falhas_vies) + " de " + str(total_m), "<= 2", falhas_vies <= 2)
    linha(verif, 7, "dispersão contra sqrt((1-A^2)/T)", "fora de 5%: " + str(falhas_razao) + " de " + str(total_m), "<= 2", falhas_razao <= 2)
    razoes = []
    ok = True
    for nome in sigmas:
        r = sigmas[nome][0] / sigmas[nome][2]
        razoes.append(round(r, 2))
        if abs(r - 10) > 1:
            ok = False
    linha(verif, 8, "escala 1/raiz(N): desvio(1e3)/desvio(1e5)", str(razoes), "perto de 10 (+-1)", ok)
    # inclinação da reta log(sigma) contra log(N) nos 4 pontos, tem que dar -1/2
    inclinacoes = []
    ok = True
    for nome in sigmas_escala:
        xs = []
        ys = []
        for j in range(len(Ns_escala)):
            xs.append(math.log(Ns_escala[j]))
            ys.append(math.log(sigmas_escala[nome][j]))
        b = inclinacao(xs, ys)
        inclinacoes.append(round(b, 4))
        if abs(b + 0.5) > 0.02:
            ok = False
    print("  inclinação de log(sigma) contra log(N) (4 pontos):", inclinacoes)
    linha(verif, 23, "inclinação de log(sigma) contra log(N) em 4 pontos (F5)", str(inclinacoes), "-0.5 +- 0.02", ok)


def s4_intervalos():
    print("S4: intervalos")
    # minhas fórmulas contra as do scipy, todo k
    maior = 0
    for n in [5, 20, 100, 1000]:
        for k in range(n + 1):
            t = stats.binomtest(k, n)
            rw = t.proportion_ci(confidence_level=0.95, method="wilson")
            rc = t.proportion_ci(confidence_level=0.95, method="exact")
            w = wilson(k, n)
            c = clopper_pearson(k, n)
            for d in [abs(w[0] - rw.low), abs(w[1] - rw.high), abs(c[0] - rc.low), abs(c[1] - rc.high)]:
                if d > maior:
                    maior = d
    print("  maior diferença contra o scipy:", maior)
    linha(verif, 9, "Wilson e Clopper-Pearson contra o scipy", "maior diferença " + str(maior), "< 1e-9", maior < 1e-9)

    # extremos
    ok = intervalo_A(0, 0, "wilson") is None
    w0 = wald(0, 10)
    if w0[0] != w0[1]:
        ok = False                      # wald com k=0 tem largura zero
    c0 = clopper_pearson(0, 10)
    if c0[0] != 0.0 or c0[1] <= 0:
        ok = False
    for met in METODOS:
        for nm, np_ in [(0, 10), (10, 0), (1, 1)]:
            A_inf, A_sup = intervalo_A(nm, np_, met)
            if A_inf < -1 or A_sup > 1:
                ok = False
    linha(verif, 10, "casos extremos", "T=0, k=0, k=n, limites em [-1, 1]", "sem divisão por zero", ok)

    # configuração de estresse: N_exp=100, A=0.8 (T ~ 10 por região)
    gs = criar_geradores(cfg["semente_mae"], 20)
    mu_m, mu_p = medias_esperadas(100, pesos, [0.8] * 10)
    acertos = [0, 0, 0]
    tent = 0
    for r in range(M):
        n_m, n_p = gerar_contagens(mu_m, mu_p, gs[15])
        for i in range(10):
            T = n_m[i] + n_p[i]
            if T == 0:
                continue
            tent = tent + 1
            for j in range(3):
                inf, sup = intervalo_p(n_m[i], T, METODOS[j])
                if inf <= 0.9 and 0.9 <= sup:
                    acertos[j] = acertos[j] + 1
    coberturas = []
    coberturas_ex = []
    min_curva_T100 = []
    for j in range(3):
        coberturas.append(round(100 * acertos[j] / tent, 1))
        coberturas_ex.append(round(100 * cobertura_media(10, 0.9, METODOS[j]), 2))
        menor_t100 = 100.0
        for q in range(101):
            A_q = -0.98 + 1.96 * q / 100
            c_q = 100 * cobertura_exata(100, (1 + A_q) / 2, METODOS[j])
            if c_q < menor_t100:
                menor_t100 = c_q
        min_curva_T100.append(round(menor_t100, 2))
    print("  estresse (T~10, A=0.8), cobertura MC:", coberturas, "| exata:", coberturas_ex, "| min curva T=100:", min_curva_T100)
    garantir_pasta("resultados")
    arq_est = open("resultados/cobertura_estresse.csv", "w", newline="")
    esc_est = csv.writer(arq_est, lineterminator="\n")
    esc_est.writerow(["metodo", "cob_estresse_mc_pct", "cob_estresse_exata_pct", "min_curva_T100_pct"])
    for j in range(3):
        esc_est.writerow([METODOS[j], coberturas[j], coberturas_ex[j], min_curva_T100[j]])
    arq_est.close()

    # F2 com Wilson e Clopper-Pearson
    mu_m, mu_p = medias_esperadas(10000, pesos, cenarios["global"])
    n_m, n_p = gerar_contagens(mu_m, mu_p, gs[16])
    A_est = []
    w_ic = []
    c_ic = []
    for i in range(10):
        A_est.append(assimetria(n_m[i], n_p[i]))
        w_ic.append(intervalo_A(n_m[i], n_p[i], "wilson"))
        c_ic.append(intervalo_A(n_m[i], n_p[i], "clopper_pearson"))
    fig_assimetrias_ic(A_est, w_ic, c_ic, cenarios["global"], "F2 - Wilson e Clopper-Pearson (global, N=10000)", "F2_assimetrias_intervalos")


def s5_cobertura():
    print("S5: cobertura (demora um pouco)")
    # A=0: cobertura exata em T=100, 1000 e 10000
    for n in [100, 1000, 10000]:
        cobs = []
        for m in METODOS:
            cobs.append(round(100 * cobertura_exata(n, 0.5, m), 2))
        print("  T =", n, cobs)
    curvas = []
    for m in METODOS:
        xs = []
        ys = []
        for j in range(101):
            A = -0.98 + 1.96 * j / 100
            xs.append(A)
            ys.append(100 * cobertura_exata(100, (1 + A) / 2, m))
        curvas.append((m, xs, ys))
    fig_cobertura_curva(100, curvas, "F6_cobertura_curva_T100")

    guardado = {}
    rotulos = []
    mc_pm = {"wald": [], "wilson": [], "clopper_pearson": []}
    ex_pm = {"wald": [], "wilson": [], "clopper_pearson": []}
    fora = 0
    contas = 0
    for nome in cenarios:
        A = cenarios[nome]
        A_tot = A_total_esperado(pesos, A)
        for N in Ns:
            todos_m, todos_p = carregar_contagens(nome_arquivo(nome, N))
            rotulos.append(nome[0:4] + " " + str(N))
            for m in METODOS:
                c_reg, tent, c_tot = cobertura_mc(todos_m, todos_p, A, pesos, m)
                e_reg = cobertura_exata_regioes(N, pesos, A, m, guardado)
                e_tot = cobertura_media(N, (1 + A_tot) / 2, m)
                erro_reg = math.sqrt(c_reg * (1 - c_reg) / tent)
                erro_tot = math.sqrt(c_tot * (1 - c_tot) / M)
                contas = contas + 2
                if abs(c_reg - e_reg) > tol["sigmas_cobertura"] * erro_reg + 0.002:
                    fora = fora + 1
                if abs(c_tot - e_tot) > tol["sigmas_cobertura"] * erro_tot + 0.002:
                    fora = fora + 1
                cob[(nome, N, m)] = (c_reg, e_reg, c_tot, e_tot)
                mc_pm[m].append(100 * c_reg)
                ex_pm[m].append(100 * e_reg)
            print("  ", nome, N, "ok")
    fig_cobertura_mc(rotulos, mc_pm, ex_pm, "F6_cobertura_mc")
    linha(verif, 11, "cobertura monte carlo contra soma exata", "fora de 3 sigma: " + str(fora) + " de " + str(contas), "<= 4", fora <= 4)


def s6_testes():
    print("S6: testes")
    for n in [100, 1000, 10000, 100000]:
        print("  falso positivo exato, n =", n, ":", round(100 * fp_exato(n, alfa), 2), "%")
    gs = criar_geradores(cfg["semente_mae"], 20)
    nulos = {}
    for j in range(3):
        nulos[Ns[j]] = gerar_nulos(Ns[j], pesos, M, gs[9 + j])

    fora = 0
    for nome in cenarios:
        A = cenarios[nome]
        A_tot = A_total_esperado(pesos, A)
        for N in Ns:
            todos_m, todos_p = carregar_contagens(nome_arquivo(nome, N))
            maxz_nulo = nulos[N][0]
            cont = [0, 0, 0, 0, 0]
            rej_reg = 0
            tent_reg = 0
            for r in range(M):
                nm = todos_m[r]
                np_ = todos_p[r]
                if pvalor_exato(soma(nm), soma(nm) + soma(np_)) <= alfa:
                    cont[0] = cont[0] + 1
                ps = pvalores_regioes(nm, np_)
                for p in ps:
                    tent_reg = tent_reg + 1
                    if p <= alfa:
                        rej_reg = rej_reg + 1
                if bonferroni_rejeita(ps, alfa):
                    cont[1] = cont[1] + 1
                if p_global(max_abs_z(nm, np_), maxz_nulo) <= alfa:
                    cont[2] = cont[2] + 1
                est0, gl0 = chi2_zero(nm, np_)
                if stats.chi2.sf(est0, gl0) <= alfa:
                    cont[3] = cont[3] + 1
                esth, glh = chi2_hom(nm, np_)
                if stats.chi2.sf(esth, glh) <= alfa:
                    cont[4] = cont[4] + 1
            taxas = []
            for x in cont:
                taxas.append(round(x / M, 4))
            rej[(nome, N)] = taxas
            # monte carlo contra o valor exato (regiões e total)
            taxa_reg = rej_reg / tent_reg
            ex_reg = 0
            for i in range(len(pesos)):
                ex_reg = ex_reg + poder_exato(int(N * pesos[i]), abs(A[i]), alfa)
            ex_reg = ex_reg / len(pesos)
            ex_tot = poder_exato(N, abs(A_tot), alfa)
            e1 = math.sqrt(max(taxa_reg * (1 - taxa_reg), 0.0005) / tent_reg)
            e2 = math.sqrt(max(taxas[0] * (1 - taxas[0]), 0.0005) / M)
            if abs(taxa_reg - ex_reg) > 3 * e1 + 0.005:
                fora = fora + 1
            if abs(taxas[0] - ex_tot) > 3 * e2 + 0.005:
                fora = fora + 1
            print("  ", nome.ljust(11), str(N).ljust(7), "total / bonf / max|z| / chi2_0 / chi2_hom:", taxas)
    linha(verif, 12, "falso positivo e poder (MC contra exato)", "fora de 3 sigma: " + str(fora) + " de 18", "<= 2", fora <= 2)

    t10 = rej[("compensacao", 10000)]
    ok = t10[0] < 0.09 and t10[4] > 0.99 and t10[3] > 0.99
    linha(verif, 13, "compensação local (N=1e4)", "total rejeita " + str(t10[0]) + ", chi2_hom rejeita " + str(t10[4]), "total ~ falso positivo; chi2_hom ~ 1", ok)

    # chi2_hom contra o chi2 de Pearson da tabela 2 x k
    todos_m, todos_p = carregar_contagens(nome_arquivo("compensacao", 10000))
    maior = 0
    for r in range(5):
        pearson = stats.chi2_contingency([todos_m[r], todos_p[r]], correction=False)[0]
        d = abs(pearson - chi2_hom(todos_m[r], todos_p[r])[0])
        if d > maior:
            maior = d
    linha(verif, 14, "chi2_hom = chi2 de Pearson 2xk", "maior diferença " + str(maior), "< 1e-8", maior < 1e-8)

    curvas = []
    for T in [100, 1000, 10000]:
        xs = []
        ex = []
        nm = []
        for j in range(41):
            A = 0.005 * j
            xs.append(A)
            ex.append(poder_exato(T, A, alfa))
            nm.append(poder_normal(A, T))
        curvas.append((str(T), xs, ex, nm))
    fig_poder(curvas, "F7_poder")


def s7_eficiencia():
    print("S7: eficiência e pesos desiguais")
    a = cfg["a_eficiencia"]
    gs = criar_geradores(cfg["semente_mae"], 20)
    ok = abs(A_obs(0, a) - a) < 1e-12 and abs(corrigir(A_obs(0.05, a), a) - 0.05) < 1e-12
    esperado = A_obs_total_esperado(10000, pesos, cenarios["global"], a)
    if abs(esperado - A_obs(0.05, a)) > 1e-12:
        ok = False
    A0 = [0.0] * 10
    resumo = []
    linhas_ef = []
    for j in range(3):
        N = Ns[j]
        mu_m, mu_p = medias_com_eficiencia(N, pesos, A0, a)
        lista = []
        corr = []
        rejeita = 0
        for r in range(M):
            n_m, n_p = gerar_contagens(mu_m, mu_p, gs[12 + j])
            A_tot = assimetria_total(n_m, n_p)
            lista.append(A_tot)
            corr.append(corrigir(A_tot, a))
            if pvalor_exato(soma(n_m), soma(n_m) + soma(n_p)) <= alfa:
                rejeita = rejeita + 1
        mm = media(lista)
        em = desvio_padrao(lista) / math.sqrt(M)
        sig = mm / sigma_A(0, N)
        print("  N =", N, "| A observado médio", round(mm, 4), "| significância aparente", round(sig, 2), "| 'descobertas' falsas", round(rejeita / M, 3), "| depois de corrigir", round(media(corr), 5))
        if abs(mm - a) > 3 * em or abs(media(corr)) > 3 * desvio_padrao(corr) / math.sqrt(M):
            ok = False
        resumo.append(round(sig, 2))
        linhas_ef.append([N, 0.0, round(mm, 5), round(sig, 2), round(rejeita / M, 4),
                          round(media(corr), 5), round(desvio_padrao(corr), 5), round(sigma_A(0, N), 5)])
    linha(verif, 15, "A_obs = (A + a)/(1 + A a), sem efeito físico e a=" + str(a), "significância aparente " + str(resumo), "viés = a antes de corrigir, zero depois", ok)
    fig_eficiencia("F8_eficiencia")

    # a conhecido só com incerteza
    N = 100000
    sigma_a = cfg["sigma_a"]
    mu_m, mu_p = medias_com_eficiencia(N, pesos, A0, a)
    lista_inc = []
    corr = []
    rejeita_inc = 0
    for r in range(M):
        n_m, n_p = gerar_contagens(mu_m, mu_p, gs[15])
        A_tot = assimetria_total(n_m, n_p)
        lista_inc.append(A_tot)
        a_med = a + sigma_a * gs[15].normal()
        corr.append(corrigir(A_tot, a_med))
        if pvalor_exato(soma(n_m), soma(n_m) + soma(n_p)) <= alfa:
            rejeita_inc = rejeita_inc + 1
    mm_inc = media(lista_inc)
    sig_inc = mm_inc / sigma_A(0, N)
    desv_inc = desvio_padrao(corr)
    prev_inc = math.sqrt(sigma_A(0, N) ** 2 + sigma_a ** 2)
    print("  a medido com", sigma_a, "de incerteza: desvio", round(desv_inc, 5), "previsto", round(prev_inc, 5))
    linhas_ef.append([N, sigma_a, round(mm_inc, 5), round(sig_inc, 2), round(rejeita_inc / M, 4), round(media(corr), 5), round(desv_inc, 5), round(prev_inc, 5)])
    garantir_pasta("resultados")
    arq_ef = open("resultados/eficiencia_resumo.csv", "w", newline="")
    esc_ef = csv.writer(arq_ef, lineterminator="\n")
    esc_ef.writerow(["N", "sigma_a", "A_obs_medio", "sig_aparente", "taxa_falsa_descoberta", "media_corrigida", "desvio_corrigido", "desvio_previsto"])
    for l_ef in linhas_ef:
        esc_ef.writerow(l_ef)
    arq_ef.close()

    # pesos desiguais: o total não é a média simples das regiões
    pd = cfg["pesos_desiguais"]
    A = cenarios["compensacao"]
    mu_m, mu_p = medias_esperadas(10000, pd, A)
    tot = []
    simp = []
    for r in range(M):
        n_m, n_p = gerar_contagens(mu_m, mu_p, gs[17])
        tot.append(assimetria_total(n_m, n_p))
        simp.append(media_simples(n_m, n_p))
    print("  pesos desiguais: A_total esperado", round(A_total_esperado(pd, A), 4), "| média do total", round(media(tot), 4), "| média simples das regiões", round(media(simp), 4))


def s8_tabelas():
    print("S8: tabela de métricas")
    garantir_pasta("resultados")
    arq = open("resultados/metricas.csv", "w", newline="")
    esc = csv.writer(arq, lineterminator="\n")
    esc.writerow(["cenario", "N", "A_total", "vies_total", "erro_vies", "desvio_total", "previsto", "razao", "pull_media", "pull_desvio", "cob_wald_reg", "cob_wilson_reg", "cob_cp_reg", "cob_wald_tot", "cob_wilson_tot", "cob_cp_tot", "rej_total", "rej_bonferroni", "rej_maxz", "rej_chi2_0", "rej_chi2_hom"])
    for nome in cenarios:
        for N in Ns:
            m = metricas[(nome, N)][len(pesos)]
            linha_csv = [nome, N, round(m["A_verdadeiro"], 4), round(m["vies"], 6), round(m["erro_vies"], 6), round(m["desvio"], 6), round(m["previsto"], 6), round(m["razao"], 4), round(m["pull_media"], 4), round(m["pull_desvio"], 4)]
            for met in METODOS:
                linha_csv.append(round(100 * cob[(nome, N, met)][0], 2))
            for met in METODOS:
                linha_csv.append(round(100 * cob[(nome, N, met)][2], 2))
            linha_csv = linha_csv + rej[(nome, N)]
            esc.writerow(linha_csv)
    arq.close()
    print("  resultados/metricas.csv")


def s9_espaco_fase():
    print("S9: espaço de fase (etapa 2)")
    ef = cfg["espaco_fase"]
    m_min = m1 + m2
    m_max = M_B - m3
    print("  m(p pbar) vai de", round(m_min, 4), "a", round(m_max, 4), "GeV")
    lo, hi = s23_limites(m_min ** 2)
    ok1 = hi - lo < 1e-6
    pw = pesos_por_faixa(ef["n_faixas"], 200)
    pw_fino = pesos_por_faixa(ef["n_faixas"], 1000)
    dif = 0
    for i in range(len(pw)):
        dif = max(dif, abs(pw[i] - pw_fino[i]))
    pesos_arredondados = []
    for w in pw:
        pesos_arredondados.append(round(w, 4))
    print("  pesos w_i:", pesos_arredondados)
    linha(verif3, 18, "pesos somam 1 e não mudam com mais passos", "soma " + str(round(soma(pw), 9)) + ", maior mudança " + str(dif), "soma = 1, mudança < 1e-4", abs(soma(pw) - 1) < 1e-9 and dif < 1e-4 and ok1)
    u, b = largura_otima()
    melhor_b = 0
    melhor_z = 0
    for j in range(1, 400):
        bb = 0.03 * j
        z = significancia_esperada(bb)
        if z > melhor_z:
            melhor_z = z
            melhor_b = bb
    print("  largura ótima: b*/s =", round(b, 3), "(varredura analítica", round(melhor_b, 2), ")")
    linha(verif3, 19, "largura ótima b*/s (bisseção e varredura)", "b*/s = " + str(round(b, 3)) + ", varredura " + str(round(melhor_b, 2)), "perto de 2.8", 2.78 < b < 2.82 and abs(melhor_b - b) < 0.05)

    # varredura por monte carlo: densidade uniforme numa faixa de 40 s, bump gaussiano no meio
    gs = criar_geradores(cfg["semente_mae"], 30)
    s = ef["s"]
    L = 40 * s
    N = 20000
    bs = []
    z_mc = []
    z_teo = []
    for b_sobre_s in [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]:
        bb = b_sobre_s * s
        T_b = N * bb / L
        A_b = ef["A0"] * s * math.sqrt(2 * math.pi) / bb * math.erf(bb / (2 * math.sqrt(2) * s))
        mu_m = T_b * (1 + A_b) / 2
        mu_p = T_b * (1 - A_b) / 2
        zs = []
        for r in range(300):
            nm = int(gs[18].poisson(mu_m))
            np_ = int(gs[18].poisson(mu_p))
            zs.append((nm - np_) / math.sqrt(nm + np_))
        bs.append(b_sobre_s)
        z_mc.append(media(zs))
        z_teo.append(significancia_esperada(b_sobre_s))
    i_max = z_mc.index(max(z_mc))
    escala = z_mc[i_max] / z_teo[i_max]
    z_teo_esc = []
    for z in z_teo:
        z_teo_esc.append(z * escala)
    fig_largura(bs, z_mc, z_teo_esc, b, "F9_largura")
    print("  melhor largura na varredura MC: b/s =", bs[i_max])

    # o efeito passando pelo pipeline da etapa 1, com os pesos do espaço de fase
    centros = centros_das_faixas(ef["n_faixas"])
    A = []
    for m in centros:
        A.append(perfil_gaussiano(m, ef["A0"], ef["m0"], ef["s"]))
    mu_m, mu_p = medias_esperadas(10000, pw, A)
    n_m, n_p = gerar_contagens(mu_m, mu_p, gs[19])
    est, gl = chi2_hom(n_m, n_p)
    print("  A_total esperado", round(A_total_esperado(pw, A), 4), "| maior A_i", round(max(A), 3),
          "| p do chi2_hom", stats.chi2.sf(est, gl))
    arq = open("resultados/etapa2_pesos.csv", "w", newline="")
    esc = csv.writer(arq, lineterminator="\n")
    esc.writerow(["faixa", "m_centro_GeV", "peso", "A_i"])
    for i in range(len(pw)):
        esc.writerow([i + 1, round(centros[i], 4), round(pw[i], 5), round(A[i], 4)])
    arq.close()


def s10_dalitz():
    print("S10: Dalitz (etapa 3)")
    dz = cfg["dalitz"]
    nx = dz["nx"]
    ny = dz["ny"]
    minimo = dz["min_eventos"]
    gs = criar_geradores(cfg["semente_mae"], 30)
    pos, ef = gerar_posicoes(dz["N_eventos"], gs[20])
    print("  fração aceita na caixa:", round(ef, 3))

    # histograma das posições contra os pesos da etapa 2
    pw = pesos_por_faixa(cfg["espaco_fase"]["n_faixas"], 200)
    m_min = m1 + m2
    larg = (M_B - m3 - m_min) / 10
    conta = [0] * 10
    for (s12, s23) in pos:
        f = int((math.sqrt(s12) - m_min) / larg)
        if f == 10:
            f = 9
        conta[f] = conta[f] + 1
    maior_z = 0
    N = len(pos)
    for f in range(10):
        z = abs(conta[f] / N - pw[f]) / math.sqrt(pw[f] * (1 - pw[f]) / N)
        if z > maior_z:
            maior_z = z
    linha(verif3, 20, "pesos do espaço de fase x histograma das posições", "maior desvio " + str(round(maior_z, 2)) + " sigma", "< 4", maior_z < 4)

    maior = 0
    for i in range(1000):
        e = erro_conservacao(pos[i][0], pos[i][1])
        if e > maior:
            maior = e
    linha(verif3, 21, "conservação de E e p, massas, s12 e s23 recalculados", "maior erro " + str(maior), "< 1e-9", maior < 1e-9)

    # sem efeito nenhum: S ~ N(0,1) e p do chi2_hom uniforme
    todos_S, ps = rodar_nulos(pos, gs[21], dz["R_nulos"], nx, ny, minimo)
    m = media(todos_S)
    d = desvio_padrao(todos_S)
    abaixo = 0
    for p in ps:
        if p < 0.05:
            abaixo = abaixo + 1
    ordenados = sorted(ps)
    dist = 0
    for i in range(len(ordenados)):
        dist = max(dist, abs((i + 1) / len(ordenados) - ordenados[i]), abs(i / len(ordenados) - ordenados[i]))
    print("  nulo: média de S", round(m, 3), "desvio", round(d, 3), "| p<0.05 em", abaixo, "de", len(ps))
    linha(verif3, 22, "nulo (A = 0): S ~ N(0,1) e p uniforme", "média " + str(round(m, 3)) + ", desvio " + str(round(d, 3)) + ", p<0.05 em " + str(abaixo) + " de " + str(len(ps)), "média ~ 0, desvio ~ 1", abs(m) < 0.05 and abs(d - 1) < 0.05 and dist < 0.2 and abaixo <= 12)
    fig_hist_S(todos_S, "F10 - S dos bins sem efeito (nulo)", "F10_hist_S_nulo")

    # com efeito: um bloco no limiar de m(p pbar) e uma versão com um segundo bloco oposto
    um_bloco = dz["um_bloco"]
    dois_blocos = dz["dois_blocos"]
    garantir_pasta("resultados")
    arq_dz = open("resultados/etapa3_dalitz.csv", "w", newline="")
    esc_dz = csv.writer(arq_dz, lineterminator="\n")
    esc_dz.writerow(["cenario", "n_bins", "A_total", "p_total", "p_chi2_hom", "menor_S", "maior_S", "maior_abs_S"])
    for nome, blocos in [("um_bloco", um_bloco), ("dois_blocos", dois_blocos)]:
        cargas = atribuir_cargas(pos, blocos, gs[22])
        cx, cy, lm, lp = contar_por_bin(pos, cargas, nx, ny, minimo)
        A_tot = (soma(lm) - soma(lp)) / (soma(lm) + soma(lp))
        p_tot = pvalor_exato(soma(lm), soma(lm) + soma(lp))
        p_hom = p_do_chi2_hom(lm, lp)
        S = S_com_alfa(lm, lp)
        maior_abs = max(max(S), -min(S))
        print("  ", nome, "| A total", round(A_tot, 4), "| p do total", round(p_tot, 6), "| p do chi2_hom", p_hom, "| menor S", round(min(S), 2), "| maior S", round(max(S), 2), "| maior |S|", round(maior_abs, 2))
        esc_dz.writerow([nome, len(S), round(A_tot, 4), round(p_tot, 6), p_hom, round(min(S), 2), round(max(S), 2), round(maior_abs, 2)])
        fig_mapa_S(cx, cy, S, "F10 - mapa de S (" + nome + ")", "F10_mapa_" + nome)
    arq_dz.close()


def s11_emv():
    print("S11: máxima verossimilhança e Fisher")
    maior_p = 0
    maior_v = 0
    for n in [20, 100, 1000]:
        for fr in [0.3, 0.5, 0.7]:
            k = int(round(fr * n))
            p_num = emv_binomial(k, n)
            maior_p = max(maior_p, abs(p_num - k / n))
            p = k / n
            v_analitica = 1 / informacao_fisher(p, n)
            v_num = variancia_numerica(k, n)
            maior_v = max(maior_v, abs(v_num / v_analitica - 1))
    print("  maior diferença no EMV:", maior_p, "| na variância 1/I:", maior_v)
    linha(verif, 17, "EMV binomial = k/n e 1/I = p(1-p)/n", "diferenças " + str(maior_p) + " e " + str(maior_v), "< 1e-6", maior_p < 1e-6 and maior_v < 1e-6)
    # quatro intervalos pro mesmo dado
    for k, n in [(50, 100), (3, 20), (0, 20)]:
        resultado = []
        for met in ["wald", "wilson", "clopper_pearson", "verossimilhanca"]:
            inf, sup = intervalo_p(k, n, met)
            resultado.append([round(inf, 4), round(sup, 4)])
        print("  k =", k, "n =", n, "->", resultado)


def numero_do_item(l):
    return l[0]


def fechar():
    garantir_pasta("resultados")
    for nome, lista in [("resultados/verificacao.csv", verif), ("resultados/verificacao_etapa3.csv", verif3)]:
        arq = open(nome, "w", newline="", encoding="utf-8")
        esc = csv.writer(arq, lineterminator="\n")
        esc.writerow(["item", "descricao", "valor", "criterio", "passou"])
        for l in sorted(lista, key=numero_do_item):
            esc.writerow(l)
        arq.close()
    reprovadas = 0
    for l in verif + verif3:
        if l[4] != "sim":
            reprovadas = reprovadas + 1
            print("  linha que não passou:", l[0], l[1])
    print("linhas da verificação que não passaram:", reprovadas)


if __name__ == "__main__":
    inicio = time.time()
    for bloco in [s1_geracao, s2_estimadores, s3_monte_carlo, s4_intervalos, s5_cobertura, s6_testes, s7_eficiencia, s8_tabelas, s9_espaco_fase, s10_dalitz, s11_emv]:
        bloco()
    fechar()
    print("terminou em", round(time.time() - inicio), "s")
