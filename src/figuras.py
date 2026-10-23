# figuras F1 a F10. cada função salva um png em figuras/
import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def salvar(nome):
    if not os.path.exists("figuras"):
        os.makedirs("figuras")
    plt.tight_layout()
    plt.savefig("figuras/" + nome + ".png", dpi=200)
    plt.close()


def fig_contagens(n_m, n_p, titulo, nome):
    # F1
    xm = []
    xp = []
    for i in range(len(n_m)):
        xm.append(i + 1 - 0.2)
        xp.append(i + 1 + 0.2)
    plt.figure(figsize=(7, 3.5))
    plt.bar(xm, n_m, width=0.4, label="B-")
    plt.bar(xp, n_p, width=0.4, label="B+")
    plt.xlabel("região")
    plt.ylabel("contagem")
    plt.title(titulo)
    plt.legend()
    salvar(nome)


def fig_assimetrias(A_est, inferior, superior, A_verd, titulo, nome):
    # F2 (barras de Wald, a primeira versão)
    xs = []
    ys = []
    e_baixo = []
    e_cima = []
    for i in range(len(A_est)):
        if A_est[i] is not None:
            xs.append(i + 1)
            ys.append(A_est[i])
            e_baixo.append(A_est[i] - inferior[i])
            e_cima.append(superior[i] - A_est[i])
    plt.figure(figsize=(7, 3.5))
    plt.errorbar(xs, ys, yerr=[e_baixo, e_cima], fmt="o", capsize=3, label="estimado")
    plt.plot(range(1, len(A_verd) + 1), A_verd, "r_", markersize=14, label="verdadeiro")
    plt.axhline(0, color="gray", linewidth=0.5)
    plt.xlabel("região")
    plt.ylabel("assimetria")
    plt.title(titulo)
    plt.legend()
    salvar(nome)


def fig_assimetrias_ic(A_est, wilson_ic, cp_ic, A_verd, titulo, nome):
    # F2 de novo, agora com Wilson e Clopper-Pearson
    plt.figure(figsize=(7.5, 3.8))
    for i in range(len(A_est)):
        if A_est[i] is None:
            continue
        x = i + 1
        plt.plot([x - 0.1, x - 0.1], [wilson_ic[i][0], wilson_ic[i][1]], "b-", linewidth=2)
        plt.plot([x + 0.1, x + 0.1], [cp_ic[i][0], cp_ic[i][1]], "g-", linewidth=2)
        plt.plot(x, A_est[i], "ko", markersize=4)
    plt.plot(range(1, len(A_verd) + 1), A_verd, "r_", markersize=16)
    plt.plot([], [], "b-", label="Wilson")
    plt.plot([], [], "g-", label="Clopper-Pearson")
    plt.plot([], [], "r_", label="verdadeiro")
    plt.axhline(0, color="gray", linewidth=0.5)
    plt.xlabel("região")
    plt.ylabel("assimetria")
    plt.title(titulo)
    plt.legend()
    salvar(nome)


def fig_hist_total(lista_A, A_verd, sigma, titulo, nome):
    # F3
    plt.figure(figsize=(6, 3.5))
    plt.hist(lista_A, bins=40, density=True, alpha=0.6, label="simulação")
    xs = []
    ys = []
    for j in range(200):
        x = A_verd - 4 * sigma + 8 * sigma * j / 199
        xs.append(x)
        ys.append(math.exp(-0.5 * ((x - A_verd) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi)))
    plt.plot(xs, ys, "r-", label="normal prevista")
    plt.xlabel("A total estimado")
    plt.title(titulo)
    plt.legend()
    salvar(nome)


def fig_compensacao(total, meta1, meta2, titulo, nome):
    # F4: o total fica em zero, as metades não
    plt.figure(figsize=(6.5, 3.6))
    plt.hist(meta1, bins=40, density=True, alpha=0.6, label="regiões 1-5")
    plt.hist(meta2, bins=40, density=True, alpha=0.6, label="regiões 6-10")
    plt.hist(total, bins=40, density=True, alpha=0.6, label="total")
    plt.xlabel("assimetria estimada")
    plt.title(titulo)
    plt.legend()
    salvar(nome)


def fig_sigma_vs_N(Ns, sigmas_por_cenario, nome):
    # F5
    plt.figure(figsize=(5.5, 4))
    for cen in sigmas_por_cenario:
        plt.loglog(Ns, sigmas_por_cenario[cen], "o-", label=cen)
    ref = []
    for N in Ns:
        ref.append(1 / math.sqrt(N))
    plt.loglog(Ns, ref, "k--", label="1/raiz(N)")
    plt.xlabel("N_exp")
    plt.ylabel("desvio de A_total")
    plt.legend()
    salvar(nome)


def fig_cobertura_curva(n, curvas, nome):
    # F6 (curva exata). curvas: lista de (método, lista de A, lista de coberturas em %)
    plt.figure(figsize=(7, 3.8))
    for met, xs, ys in curvas:
        plt.plot(xs, ys, label=met, linewidth=1)
    plt.axhline(95, color="k", linestyle="--", linewidth=0.7)
    plt.ylim(60, 100.5)
    plt.xlabel("assimetria verdadeira A")
    plt.ylabel("cobertura exata (%)")
    plt.title("cobertura exata, T = " + str(n))
    plt.legend()
    salvar(nome)


def fig_cobertura_mc(rotulos, mc, exata, nome):
    # F6 (monte carlo contra exata)
    plt.figure(figsize=(8, 4))
    metodos = ["wald", "wilson", "clopper_pearson"]
    cores = ["r", "b", "g"]
    for j in range(3):
        xs = []
        for i in range(len(rotulos)):
            xs.append(i + (j - 1) * 0.2)
        plt.plot(xs, mc[metodos[j]], cores[j] + "o", label=metodos[j] + " (MC)")
        plt.plot(xs, exata[metodos[j]], cores[j] + "_", markersize=14)
    plt.axhline(95, color="k", linestyle="--", linewidth=0.7)
    plt.xticks(range(len(rotulos)), rotulos, fontsize=7)
    plt.ylabel("cobertura (%)")
    plt.title("cobertura das regiões: monte carlo (pontos) e exata (traços)")
    plt.legend(fontsize=7)
    plt.xticks(rotation=25)
    plt.tight_layout()
    salvar(nome)


def fig_poder(curvas, nome):
    # F7. curvas: lista de (rótulo, lista de A, poder exato, poder normal)
    plt.figure(figsize=(6, 4))
    for rotulo, xs, exato, normal in curvas:
        plt.plot(xs, exato, label="T = " + rotulo + " (exato)")
        plt.plot(xs, normal, "--", linewidth=0.8, color="gray")
    plt.axhline(0.05, color="k", linestyle=":", linewidth=0.7)
    plt.xlabel("assimetria verdadeira A")
    plt.ylabel("poder do teste (5%)")
    plt.title("poder exato (linhas) e aproximação normal (tracejado)")
    plt.legend(fontsize=8)
    salvar(nome)


def fig_eficiencia(nome):
    # F8
    plt.figure(figsize=(9, 3.6))
    plt.subplot(1, 2, 1)
    for a in [0.0, 0.02, 0.05]:
        xs = []
        ys = []
        for j in range(61):
            A = -0.3 + 0.01 * j
            xs.append(A)
            ys.append((A + a) / (1 + A * a))
        plt.plot(xs, ys, label="a = " + str(a))
    plt.xlabel("A verdadeiro")
    plt.ylabel("A observado")
    plt.legend()
    plt.subplot(1, 2, 2)
    for a in [0.005, 0.01, 0.02]:
        Ns = []
        zs = []
        for j in range(30):
            N = 10 ** (2.5 + 3.5 * j / 29)
            Ns.append(N)
            zs.append(a * math.sqrt(N))
        plt.loglog(Ns, zs, label="a = " + str(a))
    plt.axhline(3, color="gray", linestyle=":")
    plt.axhline(5, color="gray", linestyle="--")
    plt.xlabel("N_exp (A verdadeiro = 0)")
    plt.ylabel("significância aparente")
    plt.legend()
    plt.suptitle("eficiência diferente")
    salvar(nome)


def fig_largura(bs, z_mc, z_teoria, b_estrela, nome):
    # F9: significância contra a largura da região (em unidades de s)
    plt.figure(figsize=(5.5, 3.6))
    plt.plot(bs, z_mc, "o", label="monte carlo")
    plt.plot(bs, z_teoria, "-", label="teoria (escalada)")
    plt.axvline(b_estrela, color="g", linestyle="--", label="b* = " + str(round(b_estrela, 2)) + " s")
    plt.xlabel("largura da região b/s")
    plt.ylabel("Z médio")
    plt.legend()
    salvar(nome)


def fig_mapa_S(cx, cy, S, titulo, nome):
    # F10 (mapa): escala até +-6.5 para não saturar o pico de +6.2 sigma
    plt.figure(figsize=(6, 4.5))
    pontos = plt.scatter(cx, cy, c=S, cmap="RdBu_r", vmin=-6.5, vmax=6.5, s=250, marker="s")
    plt.colorbar(pontos, label="S")
    plt.xlabel("m2(p pbar) [GeV2]")
    plt.ylabel("m2(pbar K) [GeV2]")
    plt.title(titulo)
    salvar(nome)


def fig_hist_S(valores, titulo, nome):
    # F10 (histograma de S no teste nulo)
    plt.figure(figsize=(5.5, 3.5))
    plt.hist(valores, bins=40, density=True, alpha=0.6, label="S dos bins")
    xs = []
    ys = []
    for j in range(200):
        x = -4 + 8 * j / 199
        xs.append(x)
        ys.append(math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi))
    plt.plot(xs, ys, "r-", label="N(0,1)")
    plt.xlabel("S")
    plt.title(titulo)
    plt.legend()
    salvar(nome)
