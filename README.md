# Assimetrias entre B⁻ e B⁺ com dados simulados

**Resumo.** Neste projeto, investigo o comportamento estatístico e computacional de medidas de assimetria de carga entre $B^-$ e $B^+$ usando simulações de Monte Carlo. O erro do estimador $\hat A$ cai com $1/\sqrt{N}$ (inclinação $-0{,}50$ em escala log-log). Em regimes com poucos eventos por região (~10 eventos e $A = 0{,}8$), mostro que o intervalo de Wald perde cobertura bruscamente (cobre apenas 63 %, contra 96 % do intervalo de Wilson). Para um sinal localizado em $m(p\bar{p})$, mostro que a janela ótima de integração tem largura $b^\ast/s \approx 2{,}80$ (em unidades da largura $s$ do sinal). No plano de Dalitz, simulo dois blocos com assimetrias opostas que se anulam na contagem total (valor-*p* binomial $p = 0{,}97$), mas que são facilmente revelados pelo mapa local de significância $S_j$ (pico de $+6{,}17$) e pelo teste qui-quadrado global ($p \approx 1{,}3 \times 10^{-20}$). Todo o estudo foi feito com dados simulados, acompanhado por uma tabela de verificação independente (`resultados/verificacao.csv`).

![](figuras/F10_mapa_dois_blocos.png)

**F10:** Mapa de significância local $S_j$ no plano de Dalitz $(s_{12}, s_{23})$ para o cenário `dois_blocos`. Em vermelho: excesso de $B^-$; em azul: excesso de $B^+$. Os dois blocos se cancelam na contagem integrada do plano, mas aparecem claramente separados no mapa local.

Na física de partículas, a violação de CP mede a diferença de comportamento entre matéria e antimatéria. Ela é uma das condições necessárias formuladas por Sakharov (1967) para explicar por que há mais matéria do que antimatéria no Universo. No Modelo Padrão, a fase da matriz CKM gera violação de CP (Thomson, 2013), mas essa quantidade é pequena demais (por muitas ordens de grandeza) para explicar sozinha a assimetria observada na cosmologia.

Este repositório reúne um estudo computacional em Python para testar estimadores, intervalos de confiança e testes de hipóteses aplicados a medidas de assimetria. Usei como motivação física a cinemática do decaimento de três corpos $B^\pm \to p\bar{p} K^\pm$, estudado pela colaboração LHCb (LHCb Collaboration, 2014). Os artigos experimentais citados servem apenas para definir as escalas cinemáticas e a física do processo (não utilizo dados reais de detectores): eu escolho os valores de assimetria no arquivo de configuração, sorteio as contagens de Poisson para $B^-$ ($N^-$) e $B^+$ ($N^+$) e verifico se os métodos numéricos recuperam corretamente os valores colocados, a partir da definição básica

$$A = \frac{N^- - N^+}{N^- + N^+}.$$

## Notação

| Símbolo | Significado |
|---|---|
| $N^-$, $N^+$ | Contagens de eventos observados de $B^-$ e $B^+$, sorteadas de distribuições de Poisson |
| $T = N^- + N^+$ | Contagem total observada no subespaço sob análise (total, região ou bin); $T_j$ denota o bin $j$ |
| $N$, $N_\mathrm{exp}$ | Número esperado de eventos na amostra total ($10^3, 10^4$ e $10^5$ na etapa 1; $20\,000$ na etapa 3) |
| $\lambda$ | Média do sorteio de Poisson para a variável $T$ (por exemplo, $\lambda = N w_i$ na região $i$) |
| $k$, $M$ | Número de regiões cinemáticas ($k = 10$) e repetições de Monte Carlo ($M = 2000$) |
| $w_i$ | Fração do espaço de fase correspondente à região $i$ ($w_i = 0{,}1$ na etapa 1; calculada via cinemática na etapa 2) |
| $A$, $A_i$, $A_\mathrm{total}$ | Assimetria física verdadeira (em geral, na região $i$ e agregada sobre todo o espaço) |
| $\hat A$, $\sigma_A$ | Estimador pontual da assimetria e seu respectivo desvio-padrão esperado |
| $\alpha = N^-/N^+$ | Razão empírica entre as contagens totais observadas de $B^-$ e $B^+$ |
| $S_j$ | Significância estatística local calculada no bin cinemático $j$ (etapa 3) |
| $a$, $\sigma_a$ | Assimetria instrumental de detecção, $a = (\varepsilon^- - \varepsilon^+)/(\varepsilon^- + \varepsilon^+)$, e sua incerteza |
| $A_\mathrm{obs}$ | Assimetria observada na presença de eficiências instrumentais desiguais |
| $m_0$, $s$, $b$ | Posição central e largura do sinal em $m(p\bar{p})$, e extensão da janela de integração cinemática |

O desvio-padrão esperado para o estimador $\hat A$ é expresso por

$$\sigma_A = \sqrt{\frac{1 - A^2}{T}}.$$

Nas comparações teóricas contra as médias das simulações, troco $T$ pela média esperada $\bar T \approx N_\mathrm{exp}$. Deste modo, para uma amostra com $N_\mathrm{exp} = 1000$ eventos e assimetria nula, o desvio previsto é $\sigma_A = 1/\sqrt{1000} \approx 0{,}0316$.

## Etapas e montagem do estudo

### Etapa 1: estimadores, intervalos e testes de hipóteses

Dividi o espaço em $k = 10$ regiões de probabilidades idênticas ($w_i = 0{,}1$) e repeti cada configuração $M = 2000$ vezes para três tamanhos amostrais ($N_\mathrm{exp} = 10^3, 10^4$ e $10^5$), sob três cenários:
* `sem_efeito`: $A_i = 0$ em todas as 10 regiões ($A_\mathrm{total} = 0$);
* `global`: $A_i = +0{,}05$ uniformemente em todo o espaço ($A_\mathrm{total} = +0{,}05$);
* `compensacao`: $A_i = +0{,}10$ nas regiões 1 a 5 e $A_i = -0{,}10$ nas regiões 6 a 10 ($A_\mathrm{total} = 0$).

Para os intervalos de confiança (Wald, Wilson e Clopper-Pearson a 95 % de confiança), converti a proporção binomial $p = N^-/T$ para a assimetria $A = 2p - 1$ (Brown, Cai e DasGupta, 2001). Para testes de hipóteses, comparei o teste binomial exato bilateral na contagem global com o teste qui-quadrado de homogeneidade entre regiões cinemáticas ($\chi^2_\mathrm{hom}$).

### Etapa 2: espaço de fase em $m(p\bar{p})$ e largura ótima da janela

Nesta etapa, os pesos $w_i$ deixam de ser iguais e passam a ser calculados pela cinemática relativística de três corpos. Integrei numericamente pela regra dos trapézios (Landau, Páez e Bordeianu, 2024) a quantidade $2m\,[s_{23}^{\max}(m^2) - s_{23}^{\min}(m^2)]$ (a partir de $ds_{12} = 2m\,dm$) sobre faixas contíguas em $m = m(p\bar{p}) \in [2m_p, M_B - m_K] = [1{,}8765; 4{,}7856]\text{ GeV}$, utilizando as massas nominais do Particle Data Group, normalizando para $\sum w_i = 1$ (valores salvos em `resultados/etapa2_pesos.csv`). Para estudar um efeito localizado perto do limiar de produção do par próton-antipróton, apliquei um perfil gaussiano

$$A(m) = A_0 \exp\left[-\frac{(m - m_0)^2}{2 s^2}\right],$$

centrado em $m_0 = 2{,}2\text{ GeV}$, com amplitude $A_0 = 0{,}20$ e largura $s = 0{,}15\text{ GeV}$.

### Etapa 3: mapa de significância local no plano de Dalitz

Gerei $20\,000$ eventos na região cinemática permitida do plano de Dalitz bidimensional $(s_{12}, s_{23})$ de $B^\pm \to p\bar{p} K^\pm$. Em cada bin bidimensional, calculei a significância local descontando a razão de contagens totais $\alpha = N^-/N^+$, seguindo o método de Bediaga et al. (2009):

$$S_j = \frac{N^-_j - \alpha N^+_j}{\sqrt{\alpha T_j}}.$$

No cenário `dois_blocos` (parametrizado em `config.json`), coloquei dois blocos gaussianos bidimensionais com sinais opostos no plano $(s_{12}, s_{23})$: um com $A_0 = +0{,}30$ centrado em $(4{,}5; 12{,}0)\text{ GeV}^2$ com dispersões $(1{,}2; 4{,}0)\text{ GeV}^2$, e outro mais largo com $A_0 = -0{,}15$ em $(15{,}0; 6{,}0)\text{ GeV}^2$ com dispersões $(3{,}0; 2{,}5)\text{ GeV}^2$. Eles foram calibrados para se cancelarem quase por completo na contagem total, de modo que a assimetria global simulada resulta compatível com zero ($A_\mathrm{total} = -0{,}0003$, com incerteza amostral $1/\sqrt{N} \approx 0{,}007$ para $N = 20\,000$).

## Principais resultados

### 1. Compensação local e resumo da etapa 1

No cenário `compensacao`, as regiões 1 a 5 possuem $A_i = +0{,}10$ e as regiões 6 a 10 possuem $A_i = -0{,}10$. Como todas as regiões têm o mesmo peso, no somatório total uma metade cancela a outra e a assimetria total esperada dá zero ($A_\mathrm{total} = 0$).

![](figuras/F4_compensacao.png)

**F4:** Histogramas da assimetria estimada $\hat A$ em $N = 10^4$ para as regiões 1 a 5 ($A = +0{,}10$), regiões 6 a 10 ($A = -0{,}10$) e contagem total (verde, centrada em zero).

Com $N = 10\,000$, a distribuição amostral de $\hat A_\mathrm{total}$ permanece centrada na origem e o teste binomial bilateral global rejeita a hipótese nula em apenas $(4{,}6 \pm 0{,}5)\%$ das repetições — ou seja, comporta-se exatamente como o falso positivo nominal de 5 %. Em contrapartida, o teste qui-quadrado de homogeneidade ($\chi^2_\mathrm{hom}$), que compara as 10 regiões entre si, detecta a assimetria em 100 % das vezes a partir de $N = 10\,000$. Esse modelo de brinquedo foi motivado pela inversão de sinal da assimetria entre regiões do espaço de fase e pelo mecanismo de compensação $\pi\pi \leftrightarrow KK$ reportado pelo LHCb (2023).

A tabela abaixo resume as 9 configurações simuladas (dados salvos em `resultados/metricas.csv`):

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

Nota metodológica: a coluna "Rejeição total" mede a taxa de falsos positivos (nominal de $\approx 5\%$) nos cenários `sem_efeito` e `compensacao` (onde $A_\mathrm{total} = 0$), e o poder de detecção no cenário `global`. A coluna "Rejeição $\chi^2_\mathrm{hom}$" quantifica falsos positivos em `sem_efeito` e `global` (assimetria homogênea entre bins), e poder estatístico em `compensacao`.

### 2. Comportamento do erro com $1/\sqrt{N}$ e calibração via pull

Pela propagação de incertezas para somas de variáveis de Poisson (Blitzstein e Hwang, 2019), o desvio-padrão amostral de $\hat A$ cai com $1/\sqrt{N}$, e o Teorema Central do Limite assegura a convergência para a distribuição normal conforme $N$ cresce (Cowan, 1998; Wasserman, 2004).

![](figuras/F5_sigma_vs_N.png)

**F5:** Desvio-padrão amostral de $\hat A_\mathrm{total}$ em função de $N_\mathrm{exp}$ em escala log-log para os três cenários; a linha tracejada representa a escala teórica $1/\sqrt{N}$.

Nas simulações cobrindo duas ordens de grandeza ($N = 10^3, 10^4$ e $10^5$, fator 100), a razão entre os desvios obtidos em $10^3$ e $10^5$ resultou em $9{,}76$ (`sem_efeito`), $10{,}17$ (`global`) e $10{,}04$ (`compensacao`), muito perto do fator 10 esperado. Incluindo a amostra com $N = 3000$, os coeficientes angulares das retas de regressão em escala log-log ficaram em $-0{,}494$, $-0{,}505$ e $-0{,}502$, consistentes com o valor esperado $-0{,}50 \pm 0{,}02$.

A variável de *pull*, definida por $(\hat A - A)/\hat\sigma_A$, teve média compatível com zero e desvio-padrão de $0{,}98$ no cenário `global` com $N = 10^4$ (ficando entre $0{,}98$ e $1{,}01$ nas nove configurações), o que confirma que o estimador não tem viés e que a incerteza está bem calibrada. A forma gaussiana foi avaliada visualmente pelo histograma empírico sobreposto à função densidade de probabilidade normal prevista (F3), sem aplicação de teste formal de hipótese de normalidade.

![](figuras/F3_hist_total_global.png)

**F3:** Histograma empírico de $\hat A_\mathrm{total}$ no cenário `global` ($N = 10^4$) comparado à curva normal teórica prevista (linha vermelha).

### 3. Cobertura empírica dos intervalos de confiança

![](figuras/F6_cobertura_curva_T100.png)

**F6:** Cobertura exata dos intervalos de confiança de Wald, Wilson e Clopper-Pearson em função de $A$ verdadeiro para $T = 100$; a linha tracejada horizontal assinala a cobertura nominal de 95 %.

Na curva de cobertura exata com $T = 100$, o intervalo de Wald perde cobertura de forma severa quando a assimetria se aproxima dos extremos físicos $\pm 1$, caindo para $63{,}3\%$ em $\lvert A\rvert = 0{,}98$. No teste de estresse em regime de poucos eventos por região ($N_\mathrm{exp} = 100$, com média $\bar T \approx 10$ e $A = 0{,}8$; salvo em `resultados/cobertura_estresse.csv`), a cobertura de Wald atinge apenas cerca de 63 % (com plena concordância entre a integração analítica e a simulação de Monte Carlo), enquanto o intervalo de Wilson preserva cobertura em torno de 96 % e o método exato de Clopper-Pearson alcança 99 %.

### 4. Eficiência instrumental diferente entre $B^-$ e $B^+$

Detectores físicos são feitos de matéria comum (prótons, nêutrons e elétrons). Por isso, mésons $K^-$ e $K^+$ interagem de maneira desigual com o material do detector: o $K^-$ (composto por quarks $s\bar{u}$) sofre maior taxa de absorção nuclear que o $K^+$ ($u\bar{s}$) porque o anti-quark $\bar{u}$ se aniquila com quarks de valência $u$ dos núcleons formando híperons ($\Lambda, \Sigma$), um efeito instrumental bem conhecido em física experimental de altas energias (LHCb Collaboration, 2014; Thomson, 2013). Disso decorre uma eficiência de detecção intrinsecamente desigual $\varepsilon^- < \varepsilon^+$, que simulei usando parâmetro instrumental $a = (\varepsilon^- - \varepsilon^+)/(\varepsilon^- + \varepsilon^+) = -0{,}02$ ($-2\%$). Como $N^\pm \propto \varepsilon^\pm (1 \pm A)$, a assimetria observada bruta $A_\mathrm{obs}$ e a sua respectiva transformação analítica de correção são:

$$A_\mathrm{obs} = \frac{A + a}{1 + a A}, \qquad A = \frac{A_\mathrm{obs} - a}{1 - a A_\mathrm{obs}}.$$

No caso em que não há assimetria física real ($A = 0$), a assimetria observada é puramente instrumental: $A_\mathrm{obs} = a$.

![](figuras/F8_eficiencia.png)

**F8:** À esquerda: assimetria observada $A_\mathrm{obs}$ em função de $A$ para $a = 0$, $0{,}02$ e $0{,}05$. À direita: significância espúria aparente em função de $N_\mathrm{exp}$ para amostras com $A = 0$, com $a = 0{,}005$, $0{,}01$ e $0{,}02$; as linhas pontilhada e tracejada assinalam os limiares de evidência ($3\sigma$) e descoberta ($5\sigma$).

Os resultados numéricos das simulações constam em `resultados/eficiencia_resumo.csv`:
* **Sem correção:** mantendo $A = 0$ e $a = -0{,}02$, a assimetria aparente medida converge para $\approx -0{,}02$. Como o erro estatístico encolhe com $1/\sqrt{N}$, a significância estatística do falso sinal salta de $0{,}66\sigma$ ($N = 10^3$) para $6{,}32\sigma$ em $N = 10^5$ (compatível com os valores previstos $\lvert a\rvert\sqrt{N} = 0{,}63\sigma, 2{,}00\sigma$ e $6{,}32\sigma$), induzindo a 100 % de falsas rejeições da hipótese nula (falso sinal a $>5\sigma$) quando a amostra é grande.
* **Correção com calibração exata de** $a$: a inversão analítica aplicada a $A_\mathrm{obs}$ recupera o valor verdadeiro com precisão numérica ($< 10^{-12}$). Nos sorteios de Monte Carlo com $a = -0{,}02$ conhecido, o viés médio após a correção zera até o nível do ruído estatístico da amostra (média observada de $3 \times 10^{-5}$ em $N = 10^5$).
* **Calibração com incerteza experimental em** $a$: quando o próprio parâmetro $a$ é medido com incerteza gaussiana $\sigma_a = 0{,}005$ ($0{,}5$ ponto percentual) em $N = 10^5$ (linha `sigma_a = 0.005` em `resultados/eficiencia_resumo.csv`), o desvio-padrão do estimador corrigido sobe de $0{,}00314$ ($0{,}31\%$, só estatístico) para $0{,}00596$ ($0{,}60\%$, incluindo incerteza sistemática), em excelente acordo com a combinação quadrática esperada $\sqrt{\sigma_A^2 + \sigma_a^2} = 0{,}00592$. Em amostras grandes, é a incerteza de calibração que domina o erro final da medida.

### 5. Espaço de fase em $m(p\bar{p})$ e largura ótima da janela

Integrando eventos em uma janela simétrica $[m_0 - b/2, m_0 + b/2]$, o número de eventos de sinal acumula-se proporcionalmente à função erro $\mathrm{erf}\left(b / (2\sqrt{2} s)\right)$, enquanto as flutuações estatísticas do total crescem com $\sqrt{b}$. A significância esperada comporta-se como:

$$Z(b) \propto \frac{\mathrm{erf}\left(\frac{b}{2\sqrt{2} s}\right)}{\sqrt{b}},$$

apresentando um ponto máximo em $b^\ast/s \approx 2{,}80$. Essa relação analítica é válida para um perfil gaussiano unidimensional sobre densidade cinemática aproximadamente uniforme, sem fundo e sem interferência entre amplitudes.

![](figuras/F9_largura.png)

**F9:** Significância estatística média $Z$ em função da largura relativa da janela de integração $b/s$. Círculos: Monte Carlo; linha contínua: modelo teórico analítico escalado; linha tracejada verde: ponto ótimo $b^\ast/s$.

### 6. Mapeamento cinemático no plano de Dalitz ($B^\pm \to p\bar{p} K^\pm$)

O mapa de significância bidimensional está exibido na figura F10, no início deste documento. Conforme registrado em `resultados/etapa3_dalitz.csv`, o bin cinemático mais expressivo do bloco positivo alcança significância local $S_j = +6{,}17$, enquanto o bin do bloco negativo atinge $S_j = -3{,}95$. Tratam-se de valores extremos locais em uma busca sobre dezenas de bins (efeito *look-elsewhere*); a significância global do decaimento é fornecida pelo teste $\chi^2_\mathrm{hom}$ sobre os 48 bins com ocupação $\ge 20$ eventos, atingindo valor-*p* nominal de $p \approx 1{,}3 \times 10^{-20}$ a partir da cauda assintótica contínua de $\chi^2$ com 47 graus de liberdade (confirmando com clareza que a assimetria local é real). Em contrapartida, na soma integrada sobre todo o plano de Dalitz, o cancelamento mútuo entre os dois blocos resulta em assimetria global nula ($A_\mathrm{total} = -0{,}0003$, estatisticamente indistinguível de zero frente a $1/\sqrt{N} \approx 0{,}007$) e valor-*p* binomial $p_\mathrm{total} \approx 0{,}97$ (registro completo em `resultados/etapa3_dalitz.csv`).

## Detalhes numéricos

* **Tabela de métricas e validação de poder (seção 1):** o desvio previsto analítico corresponde a $\sqrt{(1 - A^2)/\bar T}$, onde $\bar T \approx N_\mathrm{exp}$ denota o valor esperado total. A incerteza amostral da estimativa de viés (`erro_vies`) equivale a $\sigma_A/\sqrt{M} \approx 0{,}0007$, $0{,}0002$ e $0{,}00007$ para $N = 10^3, 10^4$ e $10^5$. A taxa de rejeição a 5 % de significância confronta o teste binomial bilateral sobre a soma total e o teste de homogeneidade $\chi^2_\mathrm{hom}$ entre regiões. Na validação de poder do item 12 da tabela de verificação (`verificacao.csv`), o poder analítico exato de referência é calculado admitindo contagem total $T$ estritamente fixa ($T = N w_i$), enquanto no experimento de Monte Carlo o valor de $T$ flutua de acordo com uma distribuição de Poisson. Para acomodar o efeito dessa flutuação de Poisson em regimes amostrais finitos, inclui-se uma margem de tolerância empírica de $+0{,}005$ (+0,5 ponto percentual) além da banda estatística de $3\sigma$; sem essa margem adicional, 3 das 18 configurações (especificamente em $N = 1000$, onde as flutuações relativas de Poisson são mais expressivas) ultrapassam a banda estrita de $3\sigma$ calculada para $T$ fixo.
* **Avaliação de cobertura em malha de Poisson (seção 3):** para manter a execução em cerca de 2 minutos, a soma binomial sobre $k$ foi restrita à janela de corte $\pm 7\sigma$ em torno da média $np$, e a média de Poisson sobre $T$ foi avaliada no intervalo $\pm 5\sqrt{\lambda}$. Para médias esperadas $\lambda \le 200$ (como na curva exata com $T = 100$ e nas regiões individuais com $N = 10^3$), somam-se todos os inteiros consecutivos; para $\lambda > 200$ (como no somatório total de $N = 10^3, 10^4$ e $10^5$, e nas regiões com $N \ge 10^4$), emprega-se uma malha discretizada representativa de $T$. A estabilidade numérica dessa amostragem foi verificada duplicando a malha de 21 para 42 nós: a variação observada na cobertura atinge no máximo $\approx 1{,}1$ ponto percentual no extremo assintótico mais desfavorável ($A = 0{,}98$ com Wald em $\lambda = 300$) e permanece tipicamente abaixo de $\approx 0{,}5$ ponto percentual na maior parte do espaço de parâmetros.
* **Extensão ótima da janela (seção 5):** o parâmetro ótimo analítico $b^\ast/s \approx 2{,}80$ é obtido resolvendo $4u e^{-u^2}/\sqrt{\pi} - \mathrm{erf}(u) = 0$, cuja raiz por bisseção fornece $u \approx 0{,}98994$ e $b^\ast/s = 2\sqrt{2}u \approx 2{,}79997$ (arredondado para $2{,}80$). O valor concorda perfeitamente com a varredura contínua de alta resolução ($2{,}79$) e com a amostragem discreta da simulação de Monte Carlo ($3{,}0$ em `F9_largura.png`).

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
