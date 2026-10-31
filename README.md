# Assimetrias entre B⁻ e B⁺ com dados simulados

## Resumo

Em decaimentos de três corpos de mésons $B^\pm$, como $B^\pm \to p\bar{p} K^\pm$, a diferença de contagens entre $B^-$ e $B^+$ não precisa ser igual em todo o espaço de fase. Ela pode mudar de valor e até inverter de sinal dependendo da região cinemática. Isso traz um problema prático importante: se uma parte do espaço de fase tiver excesso de $B^-$ e outra parte tiver excesso de $B^+$, somar todos os eventos em uma única contagem global faz um efeito cancelar o outro. O resultado integrado dá uma assimetria próxima de zero e esconde as diferenças locais. Por outro lado, quando dividimos a amostra em regiões menores para procurar esses efeitos locais, sobram menos eventos em cada divisão, o que exige cuidado com intervalos de confiança em baixa estatística, com a escolha da janela de corte cinemático e com diferenças de eficiência do detector.

Para estudar esse problema de um jeito controlado, montei simulações de Monte Carlo em Python que sorteiam contagens de Poisson para $B^-$ e $B^+$ em diferentes cenários (sem efeito, com assimetria constante e com regiões de sinais opostos). Com isso, testei como se comportam os estimadores, os intervalos de confiança, a largura da janela de seleção em $m(p\bar{p})$ e os testes de hipóteses locais e globais no plano de Dalitz $(s_{12}, s_{23})$.

**Principais conclusões e resultados quantitativos:**
* **Cancelamento global e detecção por regiões:** Quando regiões opostas têm assimetrias que se anulam na soma, o teste binomial aplicado ao total não vê diferença nenhuma em relação ao acaso. Já o teste qui-quadrado de homogeneidade ($\chi^2_\mathrm{hom}$) e o mapa de significâncias locais ($S_j$) conseguem enxergar a estrutura que foi colocada no modelo. No plano de Dalitz simulado com dois blocos de sinais opostos ($N = 20\,000$), a assimetria total fica compatível com zero ($A_\mathrm{total} = -0{,}0003$, valor-*p* binomial de $0{,}97$), mas a análise dividida em bins recupera os excessos locais (picos de $S_j = +6{,}17$ e $-3{,}95$) e rejeita a hipótese de distribuição uniforme com $p \approx 1{,}3 \times 10^{-20}$.
* **Intervalos de confiança com poucos eventos:** Ao dividir o espaço de fase em vários bins, a contagem $T$ por região fica pequena. Nesse regime, a aproximação normal de Wald falha bastante perto dos limites $\pm 1$, enquanto o método de Wilson mantém uma cobertura confiável. No teste com média de $\bar T \approx 10$ eventos por região e $A = 0{,}8$, o intervalo de Wald cobre o valor verdadeiro em apenas 63 % das vezes, contra 96 % do intervalo de Wilson e 99 % de Clopper-Pearson. Quando o número de eventos cresce, o erro do estimador $\hat A$ cai como $1/\sqrt{N}$ (inclinação de $-0{,}50$ no gráfico log-log).
* **Janela de corte em $m(p\bar{p})$ e eficiência de detecção:** Para um sinal concentrado na massa invariante $m(p\bar{p})$, a janela de integração que dá a maior razão sinal-ruído tem largura $b^\ast/s \approx 2{,}80$ (em unidades da largura $s$ do sinal), equilibrando a perda de sinal em janelas muito curtas e o excesso de flutuação estatística em janelas muito abertas. Além disso, uma diferença de eficiência não corrigida de $a = -0{,}02$ ($-2\%$) entre $B^-$ e $B^+$ gera um falso sinal que cresce com $\sqrt{N}$ (chegando a $6{,}32\sigma$ e 100 % de falsas rejeições em $N = 10^5$). A correção analítica remove esse viés, mas o erro final em amostras grandes passa a depender da incerteza sistemática $\sigma_a$ da calibração.

Todo o trabalho foi feito com dados simulados e acompanhado de uma tabela de verificação (`resultados/verificacao.csv`) que compara cada cálculo numérico com fórmulas analíticas ou funções do SciPy.

---

## Introdução e definição do problema

Na física de partículas, a violação de CP mede a diferença de comportamento entre matéria e antimatéria. Ela é uma das condições propostas por Sakharov (1967) para explicar por que o Universo observável é dominado por matéria. No Modelo Padrão, a violação de CP vem da fase complexa da matriz CKM (Thomson, 2013), mas o efeito previsto por esse mecanismo é muitas ordens de grandeza menor do que o necessário para explicar a matéria na cosmologia. Por isso, decaimentos de hádrons pesados continuam sendo estudados em detalhe.

Em decaimentos de três corpos de mésons $B^\pm$ sem charme, como $B^\pm \to p\bar{p} K^\pm$ (LHCb Collaboration, 2014, 2023), o processo passa por estados intermediários e ressonâncias que interferem entre si no espaço de fase bidimensional (o plano de Dalitz). Como as fases da amplitude variam com as massas invariantes, a assimetria entre $B^-$ e $B^+$ não precisa ter o mesmo valor em todos os pontos: uma região cinemática pode ter mais decaimentos de $B^-$ e outra região pode ter mais decaimentos de $B^+$.

Em qualquer amostra ou sub-região cinemática, a assimetria é medida a partir das contagens observadas de $B^-$ ($N^-$) e de $B^+$ ($N^+$):

$$A = \frac{N^- - N^+}{N^- + N^+}.$$

Se o espaço de fase for dividido em regiões $i$, cada uma contendo uma fração $w_i$ do total de eventos ($\sum w_i = 1$) e uma assimetria local $A_i$, a assimetria total integrada será a média ponderada $A_\mathrm{total} = \sum w_i A_i$. O problema central que abordo aqui vem direto dessa soma: se existirem regiões com $A_i > 0$ e regiões com $A_j < 0$, uma parte cancela a outra quando juntamos tudo. Assim, uma medida feita apenas sobre a contagem total pode dar $A_\mathrm{total} \approx 0$ mesmo quando existem diferenças fortes dentro do espaço de fase.

Este repositório usa simulações em Python para investigar esse problema. Os artigos experimentais citados servem como referência para as massas das partículas, os limites cinemáticos e a motivação física, mas **não uso dados reais de detectores**. Eu mesmo defino os valores de assimetria no arquivo `config.json`, sorteio as contagens de Poisson e testo se as ferramentas estatísticas conseguem recuperar o que foi colocado no modelo.

## Notação

| Símbolo | Significado |
|---|---|
| $N^-$, $N^+$ | Contagens de eventos de $B^-$ e $B^+$, sorteadas de distribuições de Poisson |
| $T = N^- + N^+$ | Contagem total observada no recorte analisado (amostra inteira, região $i$ ou bin $j$) |
| $N$, $N_\mathrm{exp}$ | Número esperado de eventos na amostra total ($10^3, 10^4$ e $10^5$ na etapa 1; $20\,000$ na etapa 3) |
| $\lambda$ | Média do sorteio de Poisson para a contagem $T$ ($\lambda = N w_i$ na região $i$) |
| $k$, $M$ | Número de regiões cinemáticas ($k = 10$) e número de repetições de Monte Carlo ($M = 2000$) |
| $w_i$ | Fração do espaço de fase da região $i$ ($w_i = 0{,}1$ na etapa 1; calculada pela cinemática na etapa 2) |
| $A$, $A_i$, $A_\mathrm{total}$ | Assimetria verdadeira colocada no modelo (geral, na região $i$ e somada no espaço todo) |
| $\hat A$, $\sigma_A$ | Estimador da assimetria e seu desvio-padrão esperado |
| $\alpha = N^-/N^+$ | Razão entre as contagens totais de $B^-$ e $B^+$ na amostra |
| $S_j$ | Significância estatística local no bin $j$ do plano de Dalitz (etapa 3) |
| $a$, $\sigma_a$ | Assimetria instrumental de detecção, $a = (\varepsilon^- - \varepsilon^+)/(\varepsilon^- + \varepsilon^+)$, e sua incerteza |
| $A_\mathrm{obs}$ | Assimetria observada antes de corrigir a diferença de eficiência |
| $m_0$, $s$, $b$ | Posição central e largura do perfil gaussiano em $m(p\bar{p})$, e largura da janela de corte |

O desvio-padrão esperado para o estimador $\hat A$ em uma contagem total $T$ é

$$\sigma_A = \sqrt{\frac{1 - A^2}{T}}.$$

Nas comparações teóricas com as médias das simulações, troco $T$ pela média esperada $\bar T \approx N_\mathrm{exp}$. Por exemplo, para $N_\mathrm{exp} = 1000$ eventos e assimetria nula, o desvio esperado é $\sigma_A = 1/\sqrt{1000} \approx 0{,}0316$.

---

## Metodologia e etapas da simulação

### Etapa 1: estimadores, intervalos de confiança e testes em regiões discretas

Dividi o espaço em $k = 10$ regiões de mesmo peso ($w_i = 0{,}1$) e repeti o sorteio $M = 2000$ vezes para três tamanhos de amostra ($N_\mathrm{exp} = 10^3, 10^4$ e $10^5$), em três cenários:
* `sem_efeito`: $A_i = 0$ em todas as 10 regiões ($A_\mathrm{total} = 0$), para checar a taxa de falsos positivos;
* `global`: $A_i = +0{,}05$ igual em todas as regiões ($A_\mathrm{total} = +0{,}05$);
* `compensacao`: $A_i = +0{,}10$ nas regiões 1 a 5 e $A_i = -0{,}10$ nas regiões 6 a 10 ($A_\mathrm{total} = 0$).

Para calcular os intervalos de confiança de 95 % (Wald, Wilson e Clopper-Pearson), achei primeiro os limites para a proporção binomial $p = N^-/T$ e depois passei para a assimetria usando $A = 2p - 1$ (Brown, Cai e DasGupta, 2001). Nos testes de hipóteses, comparei o teste binomial exato bilateral (aplicado apenas à soma total) com o teste qui-quadrado de homogeneidade ($\chi^2_\mathrm{hom}$), que compara as 10 regiões entre si.

### Etapa 2: cinemática em $m(p\bar{p})$ e janela ótima de seleção

Aqui os pesos $w_i$ deixam de ser iguais e passam a vir da cinemática relativística de três corpos. Integrei pela regra dos trapézios (Landau, Páez e Bordeianu, 2024) a largura disponível de espaço de fase $2m\,[s_{23}^{\max}(m^2) - s_{23}^{\min}(m^2)]$ (usando $ds_{12} = 2m\,dm$) em 10 faixas de massa invariante $m = m(p\bar{p}) \in [2m_p, M_B - m_K] = [1{,}8765; 4{,}7856]\text{ GeV}$, com as massas do Particle Data Group e normalizando para $\sum w_i = 1$ (`resultados/etapa2_pesos.csv`). Para simular um efeito concentrado perto do limiar de produção $p\bar{p}$ (motivado pelo acúmulo de eventos perto do limiar observado pelo LHCb; LHCb Collaboration, 2014), usei um perfil gaussiano

$$A(m) = A_0 \exp\left[-\frac{(m - m_0)^2}{2 s^2}\right],$$

centrado em $m_0 = 2{,}2\text{ GeV}$, com amplitude $A_0 = 0{,}20$ e largura $s = 0{,}15\text{ GeV}$, e varri a largura $b$ da janela de corte $[m_0 - b/2, m_0 + b/2]$.

### Etapa 3: mapeamento local no plano de Dalitz

Gerei $N_\mathrm{exp} = 20\,000$ eventos espalhados de modo uniforme na região cinemática permitida do plano de Dalitz $(s_{12}, s_{23}) = (m^2_{p\bar{p}}, m^2_{\bar{p}K})$ de $B^\pm \to p\bar{p} K^\pm$. Dividi o plano em uma grade 2D e calculei, para cada bin $j$ com pelo menos 20 eventos ($T_j \ge 20$), a significância local descontando a razão global $\alpha = N^-/N^+$, conforme o método de Miranda (Bediaga et al., 2009):

$$S_j = \frac{N^-_j - \alpha N^+_j}{\sqrt{\alpha T_j}}.$$

No cenário `dois_blocos`, coloquei artificialmente duas manchas gaussianas de sinais opostos no plano $(s_{12}, s_{23})$: uma positiva com $A_0 = +0{,}30$ centrada em $(4{,}5; 12{,}0)\text{ GeV}^2$ e larguras $(1{,}2; 4{,}0)\text{ GeV}^2$, e outra negativa mais larga com $A_0 = -0{,}15$ centrada em $(15{,}0; 6{,}0)\text{ GeV}^2$ e larguras $(3{,}0; 2{,}5)\text{ GeV}^2$. Escolhi esses valores justamente para que um bloco cancelasse o outro na contagem total do plano, testando se o mapa de $S_j$ e o teste $\chi^2_\mathrm{hom}$ conseguem recuperar os dois blocos.

---

## Principais resultados e discussão

### 1. Cancelamento na contagem global vs. detecção por homogeneidade

No cenário `compensacao`, as regiões 1 a 5 têm $A_i = +0{,}10$ e as regiões 6 a 10 têm $A_i = -0{,}10$. A figura **F4** mostra o que acontece com a assimetria estimada $\hat A$ nas 2000 simulações com $N = 10\,000$:

![](figuras/F4_compensacao.png)

**F4:** Histogramas das estimativas de $\hat A$ em $M = 2000$ repetições ($N = 10^4$) no cenário `compensacao`, mostrando separadamente o grupo das regiões 1 a 5 ($A = +0{,}10$, azul), o grupo das regiões 6 a 10 ($A = -0{,}10$, laranja) e o resultado da contagem total somada (verde).

Olhando para a figura **F4**, vê-se que as distribuições das duas metades do espaço de fase formam dois picos bem separados, um em $+0{,}10$ e outro em $-0{,}10$, enquanto o histograma verde da soma total fica estreito e cravado no zero. Esse resultado mostra na prática o problema do cancelamento: quando somamos tudo antes de testar, os excessos de sinais opostos se anulam e o teste binomial global rejeita a hipótese nula em apenas $(4{,}6 \pm 0{,}5)\%$ das vezes, o que empata com a taxa de falso positivo de 5 % esperada para quando não há sinal. Já o teste qui-quadrado de homogeneidade ($\chi^2_\mathrm{hom}$), que olha para a diferença entre as 10 regiões, identifica a variação local em 100 % das simulações a partir de $N = 10\,000$. Essa montagem simples reproduz a ideia de compensação entre regiões do espaço de fase discutida pelo LHCb (2023).

A tabela abaixo junta os resultados das 9 configurações da Etapa 1 (salvos em `resultados/metricas.csv`):

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

*Como ler as colunas de rejeição:* a coluna "Rejeição total" mede a taxa de falsos positivos (esperada em 5 %) nos cenários `sem_efeito` e `compensacao` (onde $A_\mathrm{total} = 0$), e mede o poder de detecção no cenário `global`. Já a coluna "Rejeição $\chi^2_\mathrm{hom}$" mede falsos positivos em `sem_efeito` e `global` (onde não há diferença entre as regiões), e mede o poder de detecção em `compensacao`.

### 2. Queda do erro com $1/\sqrt{N}$ e distribuição gaussiana

Pela propagação de incertezas em contagens de Poisson (Blitzstein e Hwang, 2019), o desvio-padrão de $\hat A$ deve diminuir com $1/\sqrt{N}$, e pelo Teorema Central do Limite a distribuição das estimativas deve se aproximar de uma curva normal quando $N$ aumenta (Cowan, 1998; Wasserman, 2004). As figuras **F5** e **F3** verificam essas duas propriedades nas simulações:

![](figuras/F5_sigma_vs_N.png)

**F5:** Desvio-padrão amostral de $\hat A_\mathrm{total}$ em função do número esperado de eventos $N_\mathrm{exp}$ em escala log-log para os três cenários, junto com a reta teórica proporcional a $1/\sqrt{N}$ (linha tracejada).

Na figura **F5**, vê-se que os pontos dos três cenários caem em cima da linha tracejada teórica de $10^3$ até $10^5$ eventos. Aumentando a amostra em 100 vezes (de $10^3$ para $10^5$), o erro observado cai por um fator de $9{,}76$ (`sem_efeito`), $10{,}17$ (`global`) e $10{,}04$ (`compensacao`), bem perto do fator 10 previsto pela teoria. Fazendo o ajuste linear no gráfico log-log (incluindo o ponto em $N = 3000$), as inclinações encontradas são $-0{,}494$, $-0{,}505$ e $-0{,}502$, confirmando na prática a queda com potência $-0{,}50 \pm 0{,}02$.

![](figuras/F3_hist_total_global.png)

**F3:** Histograma normalizado das 2000 estimativas de $\hat A_\mathrm{total}$ no cenário `global` ($N = 10^4$, $A = +0{,}05$) comparado com a curva normal teórica $\mathcal{N}(A, \sigma_A^2)$ (em vermelho).

O gráfico da figura **F3** mostra que o formato do histograma simulado segue de perto o sino da curva normal vermelha centrada em $+0{,}05$ com largura $\sigma_A = 0{,}010$, sem entortar para os lados. Além disso, calculando a variável *pull* $(\hat A - A)/\hat\sigma_A$ nas nove configurações, a média deu sempre próxima de zero e o desvio-padrão ficou entre $0{,}98$ e $1{,}01$, o que mostra que o estimador não tem viés e que a fórmula usada para o erro funciona bem.

### 3. Cobertura dos intervalos de confiança com poucos eventos

Quando dividimos o espaço de fase em muitas partes, cada bin fica com poucos eventos. Para ver como os intervalos de confiança se saem nessa situação, a figura **F6** traz o cálculo exato da probabilidade de cobertura para uma contagem total de $T = 100$ ao longo de vários valores de $A$:

![](figuras/F6_cobertura_curva_T100.png)

**F6:** Cobertura exata dos intervalos de Wald, Wilson e Clopper-Pearson em função da assimetria verdadeira $A$ para $T = 100$; a linha horizontal tracejada marca a meta nominal de 95 %.

O que salta aos olhos na figura **F6** é a queda da curva azul (Wald) quando $\lvert A\rvert$ chega perto das bordas $\pm 1$, caindo para 63,3 % de cobertura em $\lvert A\rvert = 0{,}98$. Enquanto isso, a curva de Wilson (laranja) fica sempre oscilando em torno da linha tracejada de 95 %, e a de Clopper-Pearson (verde) fica toda acima de 95 % por ser mais conservadora. Isso acontece porque o intervalo de Wald usa a proporção medida na própria amostra para estimar a variância, o que dá errado quando uma das contagens ($N^-$ ou $N^+$) é muito pequena. No teste de estresse com poucos eventos por região ($N_\mathrm{exp} = 100$, dando uma média de $\bar T \approx 10$ eventos por região com $A = 0{,}8$; `resultados/cobertura_estresse.csv`), a cobertura de Wald cai para 63 %, enquanto o intervalo de Wilson mantém 96 % (e Clopper-Pearson chega a 99 %), mostrando por que é melhor usar Wilson quando há poucos eventos por bin.

### 4. Diferença de eficiência do detector e erro sistemático

Os detectores são feitos de matéria comum (prótons, nêutrons e elétrons). Por isso, partículas e antipartículas interagem de forma diferente ao atravessar o detector: o méson $K^-$ ($s\bar{u}$) tem maior chance de ser absorvido pela matéria nuclear do que o $K^+$ ($u\bar{s}$), já que o antiquark $\bar{u}$ pode se aniquilar com quarks $u$ dos núcleons formando híperons ($\Lambda, \Sigma$). Do mesmo jeito, o antipróton ($\bar{p}$) tem uma seção de choque de aniquilação nuclear $\sigma(\bar{p}N)$ bem maior do que a seção de choque do próton $\sigma(pN)$ em baixas e médias energias (LHCb Collaboration, 2014; Thomson, 2013). Isso faz as eficiências de detecção serem diferentes ($\varepsilon^- \neq \varepsilon^+$), o que é descrito pelo parâmetro instrumental $a = (\varepsilon^- - \varepsilon^+)/(\varepsilon^- + \varepsilon^+)$. Como as contagens medidas são proporcionais a $\varepsilon^\pm (1 \pm A)$, a assimetria observada $A_\mathrm{obs}$ e a fórmula para corrigi-la ficam:

$$A_\mathrm{obs} = \frac{A + a}{1 + a A}, \qquad A = \frac{A_\mathrm{obs} - a}{1 - a A_\mathrm{obs}}.$$

A figura **F8** mostra como uma eficiência desigual entorta o valor medido e infla a significância aparente se não fizermos a correção:

![](figuras/F8_eficiencia.png)

**F8:** À esquerda: assimetria observada $A_\mathrm{obs}$ em função de $A$ verdadeiro para três valores ilustrativos de $a$ ($0$, $0{,}02$ e $0{,}05$). À direita: crescimento da significância aparente ($\lvert A_\mathrm{obs}\rvert/\sigma_A$) com o aumento de $N_\mathrm{exp}$ quando não há assimetria real ($A = 0$), cruzando as linhas de $3\sigma$ (pontilhada) e $5\sigma$ (tracejada). *(Nota sobre o sinal de $a$: os gráficos foram desenhados com $a > 0$ para facilitar a leitura das curvas subindo; nas simulações numéricas usei $a = -0{,}02$, mas o efeito tem exatamente o mesmo comportamento mudando apenas o sinal).*

Na esquerda da figura **F8**, vê-se que um valor de $a \neq 0$ desloca a curva inteira para longe da diagonal $A_\mathrm{obs} = A$ (quando $A = 0$, o aparelho mede $A_\mathrm{obs} = a$). Na direita, o gráfico log-log mostra o perigo disso em amostras grandes: como o erro estatístico diminui com $1/\sqrt{N}$ mas o deslocamento instrumental $a$ é fixo, a razão sinal-ruído aparente sobe com $\lvert a\rvert\sqrt{N}$ e passa facilmente de $5\sigma$.

As simulações com $A = 0$ e $a = -0{,}02$ (`resultados/eficiencia_resumo.csv`) ilustram três situações práticas:
* **Sem corrigir:** A média medida fica em $-0{,}02$ e a significância do falso sinal sobe de $0{,}66\sigma$ ($N = 10^3$) e $1{,}96\sigma$ ($N = 10^4$) para $6{,}32\sigma$ em $N = 10^5$, gerando 100 % de falsas rejeições da hipótese nula na amostra maior.
* **Corrigindo com $a$ exato:** Aplicando a equação inversa com o valor exato $a = -0{,}02$, o cálculo devolve o valor verdadeiro sem erro numérico ($< 10^{-12}$) e zera o viés das simulações (média de $3 \times 10^{-5}$ em $N = 10^5$).
* **Corrigindo com incerteza experimental em $a$:** Quando assumimos que o próprio fator de correção $a$ foi medido com uma incerteza gaussiana $\sigma_a = 0{,}005$ (0,5 ponto percentual) em $N = 10^5$, o desvio-padrão da assimetria corrigida sobe de $0{,}00314$ (erro só estatístico) para $0{,}00596$, batendo com a soma em quadratura esperada $\sqrt{\sigma_A^2 + \sigma_a^2} = 0{,}00592$. Isso mostra que, com muitos dados, é o erro sistemático de calibração que limita a medida.

### 5. Largura ótima da janela de corte em $m(p\bar{p})$

Quando o sinal de assimetria está concentrado perto de uma massa $m_0$ com largura $s$, somar os eventos dentro de uma janela $[m_0 - b/2, m_0 + b/2]$ acumula o sinal segundo a função erro $\mathrm{erf}\left(b / (2\sqrt{2} s)\right)$, enquanto a flutuação estatística do número total de eventos dentro da janela cresce com $\sqrt{b}$. Assim, a significância esperada $Z(b)$ varia como:

$$Z(b) \propto \frac{\mathrm{erf}\left(\frac{b}{2\sqrt{2} s}\right)}{\sqrt{b}}.$$

A figura **F9** compara essa curva analítica com a significância média obtida nas simulações de Monte Carlo da Etapa 2 para diferentes larguras relativas $b/s$:

![](figuras/F9_largura.png)

**F9:** Significância estatística média $Z$ em função da largura relativa da janela $b/s$ centrada em $m_0 = 2{,}2\text{ GeV}$. Círculos azuis: média das simulações de Monte Carlo; linha contínua laranja: curva teórica proporcional a $\mathrm{erf}(b/(2\sqrt{2}s))/\sqrt{b}$; linha tracejada verde: ponto ótimo analítico $b^\ast/s \approx 2{,}80$.

O gráfico da figura **F9** mostra um pico bem definido: para janelas muito curtas ($b/s < 1{,}5$), a significância cai porque o corte joga fora muitos eventos que pertencem ao sinal; depois de passar pelo máximo em $b^\ast/s \approx 2{,}80$ (e $b/s = 3{,}0$ nos pontos simulados), abrir ainda mais a janela faz a significância cair de novo, porque as bordas distantes do centro quase não têm assimetria ($A(m) \approx 0$) e só aumentam o ruído estatístico no denominador. Isso mostra por que escolher bem a janela cinemática melhora a sensibilidade antes de fazer o teste estatístico.

### 6. Mapa de significância no plano de Dalitz e recuperação de estruturas simuladas

Na Etapa 3, analisei o plano de Dalitz $(s_{12}, s_{23})$ de $B^\pm \to p\bar{p} K^\pm$ dividido em bins bidimensionais. Como olhar para dezenas de bins ao mesmo tempo aumenta a chance de achar algum bin alto por puro sorteio (efeito *look-elsewhere*), o primeiro passo foi rodar a simulação sem nenhuma assimetria colocada no modelo (`sem_efeito`), mostrada na figura **F10 (nulo)**:

![](figuras/F10_hist_S_nulo.png)

**F10 (nulo):** Histograma normalizado das significâncias locais $S_j$ dos bins do plano de Dalitz na simulação sem efeito injetado (`sem_efeito`), comparado com a curva normal padrão $\mathcal{N}(0,1)$ (em vermelho).

Na figura **F10 (nulo)**, vê-se que, quando não há assimetria programada, os valores de $S_j$ se espalham de forma simétrica em volta do zero e acompanham a gaussiana padrão $\mathcal{N}(0,1)$, ficando dentro da faixa normal de flutuação ($\lvert S_j\rvert \lesssim 3$). Isso garante que a variável $S_j$ está bem calibrada e serve de base para comparar com o cenário onde colocamos estruturas no plano.

Em seguida, a figura **F10** mostra o mapa de cores de $S_j$ obtido no cenário `dois_blocos`, onde programei duas regiões gaussianas com assimetrias de sinais opostos ($+0{,}30$ e $-0{,}15$):

![](figuras/F10_mapa_dois_blocos.png)

**F10:** Mapa de significância local $S_j$ no plano de Dalitz $(s_{12}, s_{23})$ para o cenário simulado `dois_blocos`. Os tons em vermelho indicam excesso local de $B^-$ ($S_j > 0$) e os tons em azul indicam excesso de $B^+$ ($S_j < 0$) dentro da fronteira cinemática do decaimento.

Olhando para o mapa da figura **F10**, as duas regiões programadas aparecem bem separadas: há uma mancha vermelha forte na esquerda (baixo $s_{12}$), onde o bin de maior destaque chega a $S_j = +6{,}17$, e uma mancha azul mais espalhada na parte inferior direita (alto $s_{12}$), que chega a $S_j = -3{,}95$. Os dois valores estão muito acima da faixa de $\pm 3\sigma$ vista no histograma do teste nulo.

Quando somamos todos os $20\,000$ eventos do plano de uma vez só (`resultados/etapa3_dalitz.csv`), o bloco positivo compensa o bloco negativo e a assimetria global dá $A_\mathrm{total} = -0{,}0003$ (compatível com zero dentro do erro de $1/\sqrt{N} \approx 0{,}007$), com valor-*p* binomial de $p_\mathrm{total} \approx 0{,}97$. Já quando aplicamos o teste qui-quadrado de homogeneidade aos 48 bins com $T_j \ge 20$, o teste acusa a diferença entre as regiões com valor-*p* nominal de $p \approx 1{,}3 \times 10^{-20}$ (calculado na cauda da distribuição $\chi^2$ com 47 graus de liberdade).

Dessa forma, o teste de homogeneidade identifica diferenças entre as regiões do plano de Dalitz, embora a assimetria global permaneça compatível com zero. O resultado demonstra a capacidade do modelo de detectar estruturas locais que podem ser ocultadas pela integração global, sem constituir evidência de violação de CP no processo estudado.

---

## Detalhes técnicos e limitações do modelo

### Detalhes numéricos da implementação

* **Malha de Poisson no cálculo de cobertura (Seção 3):** Para que o script completo rode em cerca de 2 minutos, cortei a soma binomial sobre $k$ em uma janela de $\pm 7\sigma$ em volta da média $np$, e avaliei a média de Poisson sobre $T$ no intervalo $\pm 5\sqrt{\lambda}$. Quando a contagem esperada é $\lambda \le 200$ (como na curva exata de $T = 100$ e nas regiões com $N = 10^3$), o código soma todos os números inteiros um por um; quando $\lambda > 200$ (somas totais e regiões com $N \ge 10^4$), usei uma malha de 21 pontos em $T$. Testei a precisão dessa malha dobrando para 42 pontos: a cobertura mudou no máximo $1{,}1$ ponto percentual no ponto mais extremo ($A = 0{,}98$ com Wald em $\lambda = 300$) e menos de $0{,}5$ ponto percentual na maior parte dos casos.
* **Margem de tolerância no poder do teste binomial (Seção 1):** No item 12 da tabela de verificação (`resultados/verificacao.csv`), o cálculo analítico usado como referência para o poder assume que o total de eventos $T$ é fixo ($T = N w_i$), enquanto no Monte Carlo o valor de $T$ varia segundo um sorteio de Poisson. Para levar em conta essa flutuação extra de $T$ em amostras menores, adicionei uma folga de $+0{,}005$ (+0,5 ponto percentual) além da banda de $3\sigma$. Sem essa folga, 3 das 18 configurações (todas em $N = 1000$, onde $T$ flutua mais em termos relativos) passam um pouco do limite calculado com $T$ fixo.
* **Cálculo da janela ótima $b^\ast/s$ (Seção 5):** Encontrei o valor analítico $b^\ast/s \approx 2{,}80$ zerando a derivada de $\mathrm{erf}(u)/\sqrt{u}$, onde $u = b/(2\sqrt{2}s)$, o que resulta na equação $4u e^{-u^2}/\sqrt{\pi} - \mathrm{erf}(u) = 0$. Resolvendo pelo método da bisseção, a raiz fica em $u \approx 0{,}98994$, dando $b^\ast/s = 2\sqrt{2}u \approx 2{,}79997$ (arredondado para $2{,}80$), valor que concorda com a busca em grade fina ($2{,}79$) e com a grade usada no Monte Carlo ($3{,}0$).

### Limitações físicas das aproximações

* **Sem eventos de fundo:** As simulações geram apenas eventos de sinal, sem incluir fundo combinatório ou outros decaimentos parecidos. Se houvesse um fundo simétrico com $B$ eventos misturado a $S$ eventos de sinal, a incerteza estatística sobre a assimetria aumentaria por um fator $\sqrt{1 + B/S}$ (Thomson, 2013).
* **Eficiência uniforme no espaço de fase e assimetria de produção:** Em colisões próton-próton ($pp$) reais, a assimetria medida junta o efeito do decaimento com a diferença na taxa de produção de $B^-$ e $B^+$ e com variações de eficiência que dependem do momento das partículas (que nos experimentos são medidas usando canais de controle como $B^\pm \to J/\psi K^\pm$; LHCb Collaboration, 2023). Na Etapa 3, assumi eficiência constante ao longo do plano de Dalitz; por isso, o fator $\alpha = N^-/N^+$ na fórmula de $S_j$ já desconta qualquer diferença global uniforme antes de olhar as regiões locais.
* **Perfis gaussianos em vez de amplitudes com interferência:** Modelei as variações locais de assimetria usando funções gaussianas simples em $m(p\bar{p})$ e em $(s_{12}, s_{23})$ sobre uma densidade plana de eventos, sem montar um modelo dinâmico completo com somas coerentes de ressonâncias (Breit-Wigner), fatores de forma e fases dependendo dos ângulos de decaimento.

---

## Como reproduzir os resultados

Estrutura de pastas e arquivos:

```text
assimetrias-b-mc/
├── src/              # Módulos com as funções numéricas, estimadores, testes e cinemática
├── config.json       # Sementes fixas, parâmetros dos cenários e configurações
├── rodar_tudo.py     # Script que roda todas as simulações e gera os CSVs e figuras
├── resultados/       # Tabelas em CSV e arquivos de contagens (inclui verificacao.csv)
├── figuras/          # Gráficos gerados pelas simulações (F1 a F10)
├── notebook.ipynb    # Caderno com a análise explicada passo a passo
└── requirements.txt  # Lista de pacotes Python utilizados
```

O arquivo `notebook.ipynb` mostra o passo a passo do projeto: análise de um sorteio único, repetição das 2000 simulações, curvas de cobertura dos intervalos, cancelamento local no cenário de compensação, efeito e correção da eficiência desigual, resultados das Etapas 2 e 3 e a tabela de verificação. Como o notebook lê os dados prontos da pasta `resultados/`, é preciso rodar o script `rodar_tudo.py` antes de executar o notebook do zero.

Para gerar novamente todos os resultados e figuras (testado em Python 3.12, levando cerca de 2 minutos com as sementes de `config.json`):

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
