# Assimetrias entre B⁻ e B⁺ com dados simulados

**Resumo.** O erro de $A$ cai com $1/\sqrt{N}$ (inclinação $-0{,}50$ em escala log-log). Com cerca de 10 eventos por região e $A = 0{,}8$, o intervalo de Wald cobre apenas 63 %, contra 96 % do intervalo de Wilson. A janela de integração ótima em $m(p\bar{p})$ possui largura $b^\ast/s \approx 2{,}8$ (em unidades da largura $s$ do sinal). No plano de Dalitz, dois blocos de sinais opostos que se cancelam na contagem total (valor-*p* binomial $p = 0{,}97$) são nitidamente identificados no mapa de significância local $S_j$ (máximo $+6{,}17$) e no teste global $\chi^2_\mathrm{hom}$ ($p \approx 1{,}3 \times 10^{-20}$). Todo o estudo é realizado com dados simulados via Monte Carlo. A tabela de verificação (`resultados/verificacao.csv`) confronta cada rotina numérica com uma referência independente (solução analítica exata, cálculo manual ou rotinas consolidadas da biblioteca SciPy).

![](figuras/F10_mapa_dois_blocos.png)

**F10:** Mapa de significância local $S_j$ no plano de Dalitz $(s_{12}, s_{23})$ para o cenário `dois_blocos`. Em vermelho: excesso de $B^-$; em azul: excesso de $B^+$. Os dois blocos se cancelam na contagem integrada do plano, mas são resolvidos localmente pelo mapa.

Na física de partículas, a violação de CP quantifica a assimetria fundamental de comportamento entre matéria e antimatéria. Ela constitui uma das três condições postuladas por Sakharov (1967) para explicar a assimetria bariônica observada no Universo. No entanto, a magnitude da violação de CP proveniente da matriz CKM no Modelo Padrão (Thomson, 2013) é ordens de grandeza inferior ao necessário para justificar a prevalência cosmológica da matéria.

Este repositório reúne um estudo computacional baseado em simulações de Monte Carlo para testar estimadores pontuais, intervalos de confiança e testes estatísticos de hipóteses aplicados à medição de assimetrias. O contexto físico adotado como motivação é a cinemática do decaimento bariônico charmless de três corpos $B^\pm \to p\bar{p} K^\pm$, investigado experimentalmente pela colaboração LHCb (LHCb Collaboration, 2014). As publicações experimentais citadas servem para fundamentar as escalas cinemáticas e a fenomenologia do decaimento (o repositório não emprega dados reais de detectores): os parâmetros de assimetria são definidos nos arquivos de configuração, as contagens de $B^-$ ($N^-$) e $B^+$ ($N^+$) são amostradas a partir de distribuições de Poisson, e as análises estatísticas verificam a recuperação dos parâmetros verdadeiros a partir da definição

$$A = \frac{N^- - N^+}{N^- + N^+}.$$

## Notação

| Símbolo | Significado |
|---|---|
| $N^-$, $N^+$ | Contagens de eventos observados de $B^-$ e $B^+$, amostradas de distribuições de Poisson |
| $T = N^- + N^+$ | Contagem total observada no subespaço sob análise (total, região ou bin); $T_j$ denota o bin $j$ |
| $N$, $N_\mathrm{exp}$ | Número esperado de eventos na amostra total ($10^3, 10^4$ e $10^5$ na etapa 1; $`20\,000`$ na etapa 3) |
| $\lambda$ | Média do processo de Poisson para a variável $T$ (por exemplo, $\lambda = N w_i$ na região $i$) |
| $k$, $M$ | Número de regiões cinemáticas ($k = 10$) e repetições de Monte Carlo ($M = 2000$) |
| $w_i$ | Fração do espaço de fase correspondente à região $i$ ($w_i = 0{,}1$ na etapa 1; calculada via cinemática na etapa 2) |
| $A$, $A_i$, $A_\mathrm{total}$ | Assimetria física verdadeira (genérica, na região $i$ e agregada sobre todo o espaço) |
| $\hat A$, $\sigma_A$ | Estimador pontual da assimetria e seu respectivo desvio-padrão esperado |
| $\alpha = N^-/N^+$ | Razão empírica entre as contagens totais observadas de $B^-$ e $B^+$ |
| $S_j$ | Significância estatística local calculada no bin cinemático $j$ (etapa 3) |
| $a$, $\sigma_a$ | Assimetria instrumental de detecção, $a = (\varepsilon^- - \varepsilon^+)/(\varepsilon^- + \varepsilon^+)$, e sua incerteza |
| $A_\mathrm{obs}$ | Assimetria observada na presença de eficiências instrumentais desiguais |
| $m_0$, $s$, $b$ | Posição central e largura da ressonância em $m(p\bar{p})$, e extensão da janela de integração cinemática |

O desvio-padrão teórico esperado para o estimador $\hat A$ é expresso por

$$\sigma_A = \sqrt{\frac{1 - A^2}{T}}.$$

Nas comparações teóricas contra médias de conjuntos de simulações, $T$ é substituído pela média esperada $\bar T \approx N_\mathrm{exp}$. Deste modo, para uma amostra total com $N_\mathrm{exp} = 1000$ eventos e assimetria nula, o desvio esperado é $\sigma_A = 1/\sqrt{1000} \approx 0{,}0316$.

## Etapas e montagem do estudo

### Etapa 1: estimadores, intervalos e testes de hipóteses

O espaço de fase é particionado em $k = 10$ regiões de probabilidades idênticas ($w_i = 0{,}1$), com $M = 2000$ pseudoexperimentos simulados para cada tamanho amostral $N_\mathrm{exp} = 10^3, 10^4$ e $10^5$, sob três regimes físicos:
* `sem_efeito`: $A_i = 0$ em todas as 10 regiões ($A_\mathrm{total} = 0$);
* `global`: $A_i = +0{,}05$ uniformemente em todo o espaço de fase ($A_\mathrm{total} = +0{,}05$);
* `compensacao`: $A_i = +0{,}10$ nas regiões 1 a 5 e $A_i = -0{,}10$ nas regiões 6 a 10 ($A_\mathrm{total} = 0$).

Avaliam-se as propriedades de cobertura a 95 % de confiança para os métodos de Wald, Wilson e Clopper-Pearson, convertidos da proporção binomial $p = N^-/T$ para o parâmetro físico $A = 2p - 1$ (Brown, Cai e DasGupta, 2001). Para testes de hipóteses, contrastam-se o teste binomial exato bilateral na contagem global e o teste qui-quadrado de homogeneidade entre regiões cinemáticas ($\chi^2_\mathrm{hom}$).

### Etapa 2: espaço de fase em $m(p\bar{p})$ e otimização da janela de seleção

Os pesos cinemáticos $w_i$ passam a ser determinados rigorosamente pela cinemática relativística de três corpos. Integra-se numericamente pela regra dos trapézios (Landau, Páez e Bordeianu, 2024) a quantidade $`2m\,[s_{23}^{\max}(m^2) - s_{23}^{\min}(m^2)]`$ (a partir de $`ds_{12} = 2m\,dm`$) sobre faixas contíguas em $m = m(p\bar{p}) \in [2m_p, M_B - m_K] = [1{,}8765; 4{,}7856]\text{ GeV}$, utilizando as massas nominais do Particle Data Group, com normalização $\sum w_i = 1$ (valores tabulados em `resultados/etapa2_pesos.csv`). Para investigar um sinal localizado próximo ao limiar de produção bariônica, implementou-se um perfil de assimetria gaussiano

$$A(m) = A_0 \exp\left[-\frac{(m - m_0)^2}{2 s^2}\right],$$

centrado em $m_0 = 2{,}2\text{ GeV}$, com amplitude $A_0 = 0{,}20$ e largura $s = 0{,}15\text{ GeV}$.

### Etapa 3: mapa de significância local no plano de Dalitz

Amostram-se $`20\,000`$ eventos na região fisicamente permitida do plano de Dalitz bidimensional $(s_{12}, s_{23})$ de $B^\pm \to p\bar{p} K^\pm$. Em cada bin bidimensional calcula-se a significância local descontando a razão de contagens totais $\alpha = N^-/N^+$,

$$S_j = \frac{N^-_j - \alpha N^+_j}{\sqrt{\alpha T_j}},$$

conforme a formulação de Bediaga et al. (2009). No cenário `dois_blocos` (parametrizado em `config.json`), distribuem-se duas fontes gaussianas no plano $(s_{12}, s_{23})$: uma com $A_0 = +0{,}30$ centrada em $(4{,}5; 12{,}0)\text{ GeV}^2$ com dispersões $(1{,}2; 4{,}0)\text{ GeV}^2$, e outra com $A_0 = -0{,}15$ em $(15{,}0; 6{,}0)\text{ GeV}^2$ com dispersões $(3{,}0; 2{,}5)\text{ GeV}^2$. Por terem sinais opostos e integrais de assimetria compensatórias, a assimetria global média teórica é quase nula ($A_\mathrm{total} \approx -0{,}003$, com erro estatístico amostral $1/\sqrt{N} \approx 0{,}007$ para $`N = 20\,000`$).

## Principais resultados

### 1. Compensação local e resumo da etapa 1

**Pergunta:** um teste estatístico aplicado apenas à contagem global é capaz de identificar sinais de assimetria com cancelamento local perfeito?

No cenário `compensacao`, as regiões 1 a 5 possuem $A_i = +0{,}10$ e as regiões 6 a 10 possuem $A_i = -0{,}10$. Devido à simetria dos pesos das regiões, as contribuições individuais cancelam-se perfeitamente no somatório global, resultando em $A_\mathrm{total} = 0$.

![](figuras/F4_compensacao.png)

**F4:** Histogramas da assimetria estimada $\hat A$ em $N = 10^4$ para as regiões 1 a 5 ($A = +0{,}10$), regiões 6 a 10 ($A = -0{,}10$) e contagem total (verde, centrada em zero).

Com $`N = 10\,000`$, a distribuição amostral de $\hat A_\mathrm{total}$ permanece centrada na origem e o teste binomial bilateral global rejeita a hipótese nula em $(4{,}6 \pm 0{,}5)$ % das repetições, concordando com o nível de falso positivo nominal de 5 %. Em contraste, o teste $\chi^2_\mathrm{hom}$, comparando as 10 regiões entre si, detecta a assimetria com 100 % de poder a partir de $`N = 10\,000`$. Esse modelo foi motivado pela inversão de sinal da assimetria em diferentes domínios cinemáticos do espaço de fase e pelo mecanismo de compensação $\pi\pi \leftrightarrow KK$ reportado pelo LHCb (2023).

A tabela abaixo sintetiza os resultados das 9 configurações de simulação (salvos em `resultados/metricas.csv`):

| Cenário | $N_\mathrm{exp}$ | $A$ verdadeiro | Viés médio | Desvio obs. | Desvio prev. | Rejeição total | Rejeição $\chi^2_\mathrm{hom}$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| sem_efeito | 1000 | 0,00 | +0,0005 | 0,0310 | 0,0316 | 4,2 % | 5,0 % |
| sem_efeito | 10000 | 0,00 | −0,0002 | 0,0100 | 0,0100 | 4,9 % | 5,3 % |
| sem_efeito | 100000 | 0,00 | 0,0000 | 0,0032 | 0,0032 | 4,8 % | 4,9 % |
| global | 1000 | 0,05 | +0,0004 | 0,0316 | 0,0316 | 34,7 % | 4,9 % |
| global | 10000 | 0,05 | 0,0000 | 0,0098 | 0,0100 | 99,95 % | 5,4 % |
| global | 100000 | 0,05 | −0,0001 | 0,0031 | 0,0032 | 100 % | 5,4 % |
| compensacao | 1000 | 0,00 | +0,0002 | 0,0318 | 0,0316 | 5,1 % | 55,4 % |
| compensacao | 10000 | 0,00 | −0,0002 | 0,0099 | 0,0100 | 4,6 % | 100 % |
| compensacao | 100000 | 0,00 | −0,0001 | 0,0032 | 0,0032 | 5,1 % | 100 % |

Nota metodológica: a coluna "Rejeição total" mede taxa de falsos positivos (nominal de ≈ 5 %) nos cenários `sem_efeito` e `compensacao` (onde $A_\mathrm{total} = 0$), e poder de detecção no cenário `global`. A coluna "Rejeição $\chi^2_\mathrm{hom}$" quantifica a taxa de falsos positivos em `sem_efeito` e `global` (assimetria uniforme entre bins), e poder estatístico em `compensacao`.

### 2. Comportamento assintótico com $1/\sqrt{N}$ e calibração via pull

Conforme a teoria de propagação de incertezas para somas de variáveis aleatórias de Poisson (Blitzstein e Hwang, 2019), o desvio-padrão amostral de $\hat A$ escala com $1/\sqrt{N}$, e o Teorema Central do Limite garante a convergência assintótica para a distribuição normal (Cowan, 1998; Wasserman, 2004).

![](figuras/F5_sigma_vs_N.png)

**F5:** Desvio-padrão amostral de $\hat A_\mathrm{total}$ em função de $N_\mathrm{exp}$ em escala log-log para os três cenários; a linha tracejada representa a escala teórica $1/\sqrt{N}$.

Nas simulações principais cobrindo duas ordens de grandeza ($N = 10^3, 10^4$ e $10^5$, fator 100), a razão entre os desvios amostrais obtidos em $10^3$ e $10^5$ resultou em $9{,}76$ (`sem_efeito`), $10{,}17$ (`global`) e $10{,}04$ (`compensacao`), em excelente concordância com o fator teórico 10. Incluindo o ponto intermediário $N = 3000$, os coeficientes angulares das retas de regressão em escala log-log resultaram em $-0{,}494$, $-0{,}505$ e $-0{,}502$, consistentes com o valor analítico $-0{,}50 \pm 0{,}02$.

A distribuição de *pull*, definida por $(\hat A - A)/\hat\sigma_A$, obteve média compatível com zero e desvio-padrão de $0{,}98$ no cenário `global` com $N = 10^4$ (oscilando entre $0{,}98$ e $1{,}01$ nas nove configurações), comprovando ausência de viés sistemático e calibração fidedigna da incerteza estimada. O perfil gaussiano assintótico foi confirmado visualmente pelo histograma empírico sobreposto à função densidade de probabilidade normal prevista (F3), sem a imposição de teste formal de hipótese de normalidade.

![](figuras/F3_hist_total_global.png)

**F3:** Histograma empírico de $\hat A_\mathrm{total}$ no cenário `global` ($N = 10^4$) comparado à curva analítica gaussiana prevista (linha vermelha).

### 3. Cobertura empírica de intervalos de confiança

![](figuras/F6_cobertura_curva_T100.png)

**F6:** Cobertura exata dos intervalos de confiança de Wald, Wilson e Clopper-Pearson em função de $A$ verdadeiro para $T = 100$; a linha tracejada horizontal assinala a cobertura nominal de 95 %.

Conforme ilustrado na curva de cobertura exata com $T = 100$, o intervalo de Wald sofre degradação severa nas proximidades dos limites físicos $A \to \pm 1$, decaindo para 63,3 % de cobertura em $\lvert A\rvert = 0{,}98$. No teste de estresse em regime de baixa estatística por região ($N_\mathrm{exp} = 100$, média $\bar T \approx 10$ e $A = 0{,}8$; arquivado em `resultados/cobertura_estresse.csv`), a cobertura de Wald atinge apenas cerca de 63 % (concordando plenamente entre a integração analítica exata e a simulação de Monte Carlo), ao passo que o intervalo de Wilson preserva cobertura em cerca de 96 % e o método exato de Clopper-Pearson alcança 99 %.

### 4. Assimetria de eficiência instrumental entre $B^-$ e $B^+$

**Pergunta:** uma assimetria espúria induzida por detecção desigual de partículas pode ser matematicamente corrigida?

Detectores físicos são construídos de matéria barionica ordinária. Consequentemente, mésons $K^-$ e $K^+$ interagem diferentemente com a matéria nuclear (o $K^-$ contendo um quark $\bar{u}$ sofre maior taxa de absorção por formação de hiperons; Thomson, 2013). Disso decorre uma eficiência intrinsecamente desigual $\varepsilon^- < \varepsilon^+$, resultando em parâmetro instrumental $a = (\varepsilon^- - \varepsilon^+)/(\varepsilon^- + \varepsilon^+) < 0$. Adotou-se $a = -0{,}02$ (−2 %) na modelagem. Uma vez que $N^\pm \propto \varepsilon^\pm (1 \pm A)$, a assimetria observada bruta $A_\mathrm{obs}$ e a sua respectiva transformação analítica de correção são expressas por

$$A_\mathrm{obs} = \frac{A + a}{1 + a A}, \qquad A = \frac{A_\mathrm{obs} - a}{1 - a A_\mathrm{obs}}.$$

No caso de ausência de assimetria física ($A = 0$), a assimetria observada reduz-se exatamente a $A_\mathrm{obs} = a$.

![](figuras/F8_eficiencia.png)

**F8:** À esquerda: assimetria observada $A_\mathrm{obs}$ em função de $A$ para $a = 0$, $0{,}02$ e $0{,}05$. À direita: significância espúria aparente em função de $N_\mathrm{exp}$ para amostras com $A = 0$, com $a = 0{,}005$, $0{,}01$ e $0{,}02$; as linhas pontilhada e tracejada assinalam os limiares de evidência ($3\sigma$) e descoberta ($5\sigma$).

Os resultados numéricos das simulações constam em `resultados/eficiencia_resumo.csv`:
* **Sem correção:** mantendo $A = 0$ e $a = -0{,}02$, a assimetria aparente medida converge para $\approx -0{,}02$. Com o estreitamento do erro estatístico proporcional a $1/\sqrt{N}$, a significância estatística do falso sinal salta de $0{,}66\sigma$ ($N = 10^3$) e $1{,}96\sigma$ ($N = 10^4$) para $6{,}32\sigma$ em $N = 10^5$ (compatível com os valores previstos $\lvert a\rvert\sqrt{N} = 0{,}63\sigma, 2{,}00\sigma$ e $6{,}32\sigma$), produzindo 100 % de falsas descobertas em $N = 10^5$.
* **Correção com calibração exata de** $a$: a inversão analítica aplicada a $A_\mathrm{obs}$ recupera o valor verdadeiro com precisão de ponto flutuante ($< 10^{-12}$). Nos sorteios de Monte Carlo com $a = -0{,}02$ calibrado perfeitamente, o viés médio pós-correção anula-se até o ruído estocástico da amostra (média observada de $3 \times 10^{-5}$ em $N = 10^5$).
* **Calibração com incerteza experimental em** $a$: quando o parâmetro instrumental é determinado com incerteza gaussiana de calibração $\sigma_a = 0{,}005$ ($0{,}5$ ponto percentual) em $N = 10^5$ (linha `sigma_a = 0.005` em `resultados/eficiencia_resumo.csv`), o desvio-padrão do estimador corrigido eleva-se de $0{,}00314$ (0,31 %, puramente estatístico) para $0{,}00596$ (0,60 %, incluindo incerteza sistemática), em perfeita consonância com a combinação quadrática esperada $\sqrt{\sigma_A^2 + \sigma_a^2} = 0{,}00592$. Em grandes conjuntos de dados, o erro experimental torna-se dominado pela calibração sistemática.

### 5. Distribuição cinemática em $m(p\bar{p})$ e largura ótima da janela

**Pergunta:** qual extensão de janela de massa invariante em torno de um sinal gaussiano maximiza a significância estatística da descoberta?

![](figuras/F9_largura.png)

**F9:** Significância estatística média $Z$ em função da largura relativa da janela de integração $b/s$. Círculos: Monte Carlo; linha contínua: modelo teórico analítico escalado; linha tracejada verde: ponto ótimo $b^\ast/s$.

Integrando eventos em um intervalo simétrico $[m_0 - b/2, m_0 + b/2]$, o número de eventos de sinal acumula-se proporcionalmente à função erro $\mathrm{erf}\left(b / (2\sqrt{2} s)\right)$, ao passo que as flutuações estatísticas do fundo ou contagem total crescem com $\sqrt{b}$. A significância estatística esperada comporta-se como

$$Z(b) \propto \frac{\mathrm{erf}\left(\frac{b}{2\sqrt{2} s}\right)}{\sqrt{b}},$$

apresentando um máximo global bem definido em $b^\ast/s \approx 2{,}8$. Essa relação analítica é válida para perfis gaussianos unidimensionais imersos em densidade cinemática aproximadamente uniforme, na ausência de interferências de amplitudes.

### 6. Mapeamento cinemático no plano de Dalitz ($B^\pm \to p\bar{p} K^\pm$)

**Pergunta:** o cálculo de significâncias locais $S_j$ consegue isolar regiões com violação de CP quando o sinal integrado é nulo?

O mapa de significância bidimensional está exibido na figura F10, no início deste documento. De acordo com as métricas registradas em `resultados/etapa3_dalitz.csv`, o bin cinemático mais expressivo do bloco positivo alcança significância local $S_j = +6{,}17$, enquanto o bin do bloco negativo atinge $S_j = -3{,}95$. Tratam-se de valores extremos locais sujeitos ao efeito de múltiplas comparações (*look-elsewhere effect*); a significância global do decaimento é fornecida pelo teste $\chi^2_\mathrm{hom}$ sobre os 48 bins com ocupação $\ge 20$ eventos, atingindo valor-*p* de $p \approx 1{,}3 \times 10^{-20}$. Em contrapartida, na soma agregada de todo o plano, o cancelamento mútuo resulta em assimetria global nula ($A_\mathrm{total} = -0{,}0003$, perfeitamente compatível com a previsão $-0{,}003 \pm 0{,}007$) e valor-*p* binomial $p_\mathrm{total} \approx 0{,}97$ (valor exato arquivado em `resultados/etapa3_dalitz.csv`).

## Detalhes numéricos

* **Tabela de métricas (seção 1):** o desvio previsto analítico corresponde a $\sqrt{(1 - A^2)/\bar T}$, onde $\bar T \approx N_\mathrm{exp}$ denota o valor esperado total. A incerteza amostral da estimativa de viés (`erro_vies`) equivale a $\sigma_A/\sqrt{M} \approx 0{,}0007$, $0{,}0002$ e $0{,}00007$ para $N = 10^3, 10^4$ e $10^5$. A taxa de rejeição a 5 % de significância confronta o teste binomial bilateral sobre a soma total (com valor-*p* definido pelo dobro da cauda acumulada de $\min(N^-, N^+)$, truncado em 1) e o teste de homogeneidade $\chi^2_\mathrm{hom}$ entre regiões.
* **Avaliação de cobertura em malha de Poisson (seção 3):** com o intuito de viabilizar a execução numérica completa em cerca de 2 minutos, a soma binomial sobre $k$ foi restrita a uma janela de corte $\pm 7\sigma$ em torno da média $np$, e a média de Poisson sobre $T$ foi avaliada no intervalo $\pm 5\sqrt{\lambda}$. Para $\lambda \le 200$, somam-se todos os inteiros consecutivos; para $\lambda > 200$, emprega-se uma malha discretizada representativa de $T$. A integração analítica exata é empregada integralmente para $T = 100$ e $T \approx 10$; a amostragem em malha atua unicamente na validação assintótica em $N = 10^4$ e $10^5$ (`F6_cobertura_mc.png` e item 11 da tabela de verificação, onde 0 de 54 pontos divergem além de $3\sigma$). Dobrar a densidade da malha (de 21 para 42 nós) altera o valor da cobertura em no máximo $1{,}1$ ponto percentual no caso mais extremo ($A = 0{,}98$ com Wald em $\lambda = 300$) e em menos de $0{,}35$ ponto percentual nas demais configurações.
* **Extensão ótima da janela (seção 5):** o parâmetro ótimo $b^\ast/s \approx 2{,}80$ foi verificado independentemente por três abordagens numéricas concordantes: $2{,}798$ via algoritmo de bisseção da derivada analítica, $2{,}79$ via varredura contínua de alta resolução e $3{,}0$ na grade discreta da simulação de Monte Carlo (`F9_largura.png`).

## Limitações do modelo

* **Eventos de fundo:** a simulação não incorpora componentes de fundo combinatório ou de reflexão física. Para um fundo estatisticamente simétrico ($B$ eventos de fundo e $S$ de sinal), a incerteza estatística sobre a assimetria do sinal sofre uma amplificação pelo fator multiplicativo $\sqrt{1 + B/S}$ (Thomson, 2013).
* **Eficiência homogênea e assimetria de produção:** em colisões próton-próton ($pp$) em detectores reais, a assimetria bruta detectada combina violação direta de CP com assimetrias intrínsecas de detecção e de produção entre hádrons $B^-$ e $B^+$ (experimentalmente controladas por canais de calibração livres de assimetria de CP, tais como $B^\pm \to J/\psi K^\pm$; LHCb Collaboration, 2023). No modelo idealizado, assume-se eficiência uniforme no espaço de fase; na etapa 3, a normalização global $\alpha = N^-/N^+$ na definição de $S_j$ desconta assimetrias globais homogêneas antes de avaliar variações locais.
* **Ausência de amplitudes coerentes de interferência:** a variação espacial de assimetria foi parametrizada por perfis gaussianos fenomenológicos (unidimensionais em $m(p\bar{p})$ na etapa 2 e bidimensionais em $(s_{12}, s_{23})$ na etapa 3), sem a montagem formal de matrizes de transição com interferência entre ondas ressonantes parciais e fases fortes/fracas dependentes de momento.

## Como reproduzir os resultados

Estrutura do repositório:

```text
assimetrias-b-mc/
├── src/              # Módulos com rotinas numéricas, estimadores, testes e cinemática
├── config.json       # Sementes estocásticas fixadas, cenários e configurações físicas
├── rodar_tudo.py     # Pipeline automatizado de execução e validação da suíte de testes
├── resultados/       # Tabelas consolidadas em CSV (inclui a suíte verificacao.csv)
├── figuras/          # Gráficos (F1 a F10)
├── notebook.ipynb    # Caderno com análise passo a passo
└── requirements.txt  # Especificação de dependências Python
```

O caderno interativo `notebook.ipynb` reproduz os passos analíticos: exame de um pseudoexperimento individual, repetição das 2000 simulações de Monte Carlo, análise de cobertura de intervalos, cancelamento local por compensação, correção instrumental de eficiência, modelagem cinemática nas etapas 2 e 3 e exibição da tabela de verificação. O notebook consome os dados pré-processados na pasta `resultados/`, sendo necessária a execução prévia do script `rodar_tudo.py`.

Para executar a reprodução integral dos experimentos numéricos, tabelas e gráficos (testado em ambiente Python 3.12, tempo de execução estimado de $\approx 2$ minutos com as sementes fixadas em `config.json`):

```bash
pip install -r requirements.txt
python rodar_tudo.py
jupyter notebook notebook.ipynb
```

## Referências

* A. D. Sakharov, *Violation of CP Invariance, C Asymmetry, and Baryon Asymmetry of the Universe*, Pisma Zh. Eksp. Teor. Fiz. 5 (1967) 32 [JETP Lett. 5 (1967) 24].
* M. Thomson, *Modern Particle Physics*, Cambridge University Press, 2013.
* G. Cowan, *Statistical Data Analysis*, Oxford University Press, 1998.
* L. Wasserman, *All of Statistics: A Concise Course in Statistical Inference*, Springer, 2004.
* J. K. Blitzstein, J. Hwang, *Introduction to Probability*, 2nd ed., CRC Press, 2019.
* L. D. Brown, T. T. Cai, A. DasGupta, *Interval Estimation for a Binomial Proportion*, Statistical Science 16 (2001) 101.
* R. H. Landau, M. J. Páez, C. C. Bordeianu, *Computational Physics: Problem Solving with Python*, 4th ed., Wiley-VCH, 2024.
* I. Bediaga et al., *On a CP anisotropy measurement in the Dalitz plot*, Phys. Rev. D 80 (2009) 096006.
* LHCb Collaboration, *Direct CP violation in charmless three-body decays of B± mesons*, Phys. Rev. D 108 (2023) 012008.
* LHCb Collaboration, *Evidence for CP Violation in B+ → p p̄ K+ Decays*, Phys. Rev. Lett. 113 (2014) 141801.
