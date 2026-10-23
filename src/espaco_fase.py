# espaço de fase de B+ -> p pbar K+ (etapa 2)
import math

# Massas em GeV (valores do PDG)
M_B = 5.2793
M_P = 0.938272
M_K = 0.493677
m1 = M_P      # p
m2 = M_P      # pbar
m3 = M_K      # K+


def s23_limites(s12):
    raiz = math.sqrt(s12)
    E2 = (s12 - m1 ** 2 + m2 ** 2) / (2 * raiz)
    E3 = (M_B ** 2 - s12 - m3 ** 2) / (2 * raiz)
    p2 = math.sqrt(max(E2 ** 2 - m2 ** 2, 0.0))
    p3 = math.sqrt(max(E3 ** 2 - m3 ** 2, 0.0))
    s_min = (E2 + E3) ** 2 - (p2 + p3) ** 2
    s_max = (E2 + E3) ** 2 - (p2 - p3) ** 2
    return s_min, s_max


def pesos_por_faixa(n_faixas=10, passos=200):
    # integral de 2 m (s23_max - s23_min) dm em cada faixa, regra do trapézio
    m_min = m1 + m2
    m_max = M_B - m3
    larg = (m_max - m_min) / n_faixas
    areas = []
    for f in range(n_faixas):
        a = m_min + f * larg
        h = larg / passos
        tot = 0.0
        for j in range(passos + 1):
            m = a + j * h
            lo, hi = s23_limites(m * m)
            v = 2 * m * (hi - lo)
            if j == 0 or j == passos:
                tot = tot + v / 2
            else:
                tot = tot + v
        areas.append(tot * h)
    s = 0.0
    for x in areas:
        s = s + x
    pesos = []
    for x in areas:
        pesos.append(x / s)
    return pesos


def centros_das_faixas(n_faixas=10):
    m_min = m1 + m2
    larg = (M_B - m3 - m_min) / n_faixas
    centros = []
    for f in range(n_faixas):
        centros.append(m_min + (f + 0.5) * larg)
    return centros


def perfil_gaussiano(m, A0, m0, s):
    return A0 * math.exp(-(m - m0) ** 2 / (2 * s * s))


def g_largura(u):
    return 4 * u * math.exp(-u * u) / math.sqrt(math.pi) - math.erf(u)


def largura_otima():
    # raiz de 4u exp(-u^2)/raiz(pi) - erf(u) por bisseção; b* = 2 raiz(2) u* (em unidades de s)
    baixo = 0.1
    alto = 3.0
    for passo in range(60):
        meio = (baixo + alto) / 2
        if g_largura(baixo) * g_largura(meio) <= 0:
            alto = meio
        else:
            baixo = meio
    u = (baixo + alto) / 2
    return u, 2 * math.sqrt(2) * u


def significancia_esperada(b_sobre_s):
    # Z ~ erf(b / (2 raiz(2) s)) / raiz(b), sem as constantes
    return math.erf(b_sobre_s / (2 * math.sqrt(2))) / math.sqrt(b_sobre_s)
