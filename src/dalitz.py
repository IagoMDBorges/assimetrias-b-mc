# plano de Dalitz de B+ -> p pbar K+: eventos, binning e S por bin (etapa 3)
import math
from scipy import stats
from src.espaco_fase import M_B, m1, m2, m3, s23_limites
from src.estimadores import soma
from src.testes import chi2_hom

S12_MIN = (m1 + m2) ** 2
S12_MAX = (M_B - m3) ** 2
S23_MIN = (m2 + m3) ** 2
S23_MAX = (M_B - m1) ** 2


def gerar_posicoes(N, rng):
    # aceitação-rejeição na caixa; fica só o que está dentro da fronteira
    pos = []
    tentativas = 0
    while len(pos) < N:
        s12 = S12_MIN + (S12_MAX - S12_MIN) * rng.random()
        s23 = S23_MIN + (S23_MAX - S23_MIN) * rng.random()
        tentativas = tentativas + 1
        lo, hi = s23_limites(s12)
        if lo <= s23 and s23 <= hi:
            pos.append((s12, s23))
    return pos, N / tentativas


def quadrivetores(s12, s23):
    # no repouso do B: p1 ao longo de z, p2 no plano xz, p3 sai da conservação
    s13 = M_B ** 2 + m1 ** 2 + m2 ** 2 + m3 ** 2 - s12 - s23
    E1 = (M_B ** 2 + m1 ** 2 - s23) / (2 * M_B)
    E2 = (M_B ** 2 + m2 ** 2 - s13) / (2 * M_B)
    p1 = math.sqrt(max(E1 ** 2 - m1 ** 2, 0.0))
    p2 = math.sqrt(max(E2 ** 2 - m2 ** 2, 0.0))
    c12 = (E1 * E2 - (s12 - m1 ** 2 - m2 ** 2) / 2) / (p1 * p2)
    if c12 > 1:
        c12 = 1.0
    if c12 < -1:
        c12 = -1.0
    s12_ang = math.sqrt(1 - c12 ** 2)
    P1 = (0.0, 0.0, p1)
    P2 = (p2 * s12_ang, 0.0, p2 * c12)
    P3 = (-(P1[0] + P2[0]), -(P1[1] + P2[1]), -(P1[2] + P2[2]))
    E3 = math.sqrt(m3 ** 2 + P3[0] ** 2 + P3[1] ** 2 + P3[2] ** 2)
    return [(E1, P1), (E2, P2), (E3, P3)]


def erro_conservacao(s12, s23):
    # maior erro entre: energia total, momento total, massas e s12/s23 recalculados
    part = quadrivetores(s12, s23)
    massas = [m1, m2, m3]
    soma_E = 0.0
    soma_P = [0.0, 0.0, 0.0]
    erro_massa = 0.0
    for i in range(3):
        E = part[i][0]
        P = part[i][1]
        soma_E = soma_E + E
        for j in range(3):
            soma_P[j] = soma_P[j] + P[j]
        m_calc = E ** 2 - (P[0] ** 2 + P[1] ** 2 + P[2] ** 2) - massas[i] ** 2
        if abs(m_calc) > erro_massa:
            erro_massa = abs(m_calc)
    erros = [abs(soma_E - M_B), abs(soma_P[0]), abs(soma_P[1]), abs(soma_P[2]), erro_massa]
    E12 = part[0][0] + part[1][0]
    P12 = []
    for j in range(3):
        P12.append(part[0][1][j] + part[1][1][j])
    erros.append(abs(E12 ** 2 - (P12[0] ** 2 + P12[1] ** 2 + P12[2] ** 2) - s12))
    E23 = part[1][0] + part[2][0]
    P23 = []
    for j in range(3):
        P23.append(part[1][1][j] + part[2][1][j])
    erros.append(abs(E23 ** 2 - (P23[0] ** 2 + P23[1] ** 2 + P23[2] ** 2) - s23))
    return max(erros)


def assimetria_local(s12, s23, blocos):
    # soma de gaussianas 2D, cada bloco é (A0, s0, t0, sigma_s, sigma_t)
    t = 0.0
    for (A0, s0, t0, ss, st) in blocos:
        t = t + A0 * math.exp(-((s12 - s0) ** 2 / (2 * ss ** 2) + (s23 - t0) ** 2 / (2 * st ** 2)))
    if t > 0.95:
        t = 0.95
    if t < -0.95:
        t = -0.95
    return t


def atribuir_cargas(pos, blocos, rng):
    # B- com probabilidade (1 + A)/2. devolve -1 (B-) ou +1 (B+) pra cada evento
    cargas = []
    for (s12, s23) in pos:
        A = assimetria_local(s12, s23, blocos)
        if rng.random() < (1 + A) / 2:
            cargas.append(-1)
        else:
            cargas.append(+1)
    return cargas


def contar_por_bin(pos, cargas, nx, ny, minimo):
    # grade nx x ny na caixa; só entram os bins com pelo menos 'minimo' eventos
    dx = (S12_MAX - S12_MIN) / nx
    dy = (S23_MAX - S23_MIN) / ny
    cont = {}
    for i in range(len(pos)):
        ix = int((pos[i][0] - S12_MIN) / dx)
        iy = int((pos[i][1] - S23_MIN) / dy)
        if ix == nx:
            ix = nx - 1
        if iy == ny:
            iy = ny - 1
        if (ix, iy) not in cont:
            cont[(ix, iy)] = [0, 0]
        if cargas[i] == -1:
            cont[(ix, iy)][0] = cont[(ix, iy)][0] + 1
        else:
            cont[(ix, iy)][1] = cont[(ix, iy)][1] + 1
    cx = []
    cy = []
    lm = []
    lp = []
    for chave in cont:
        if cont[chave][0] + cont[chave][1] >= minimo:
            cx.append(S12_MIN + (chave[0] + 0.5) * dx)
            cy.append(S23_MIN + (chave[1] + 0.5) * dy)
            lm.append(cont[chave][0])
            lp.append(cont[chave][1])
    return cx, cy, lm, lp


def S_com_alfa(lista_m, lista_p):
    # tira a assimetria global antes de olhar a variação local
    alfa = soma(lista_m) / soma(lista_p)
    vals = []
    for i in range(len(lista_m)):
        T = lista_m[i] + lista_p[i]
        vals.append((lista_m[i] - alfa * lista_p[i]) / math.sqrt(alfa * T))
    return vals


def p_do_chi2_hom(lista_m, lista_p):
    est, gl = chi2_hom(lista_m, lista_p)
    return float(stats.chi2.sf(est, gl))


def rodar_nulos(pos, rng, R, nx, ny, minimo):
    # R experimentos com A = 0 no plano todo; devolve todos os S e os p do chi2_hom
    todos_S = []
    ps = []
    for e in range(R):
        cargas = atribuir_cargas(pos, [], rng)
        cx, cy, lm, lp = contar_por_bin(pos, cargas, nx, ny, minimo)
        todos_S = todos_S + S_com_alfa(lm, lp)
        ps.append(p_do_chi2_hom(lm, lp))
    return todos_S, ps
