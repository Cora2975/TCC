"""Gera Mmax por evento de Tmax diária acima do percentil 95.

O cenário de Mmax é Day-Outdoors / Young_adult / 3 h, o mesmo cenário de
destaque em ``comparacao_historico_modelo.ipynb``. Um evento é uma sequência
de pelo menos três dias de calendário com Tmax diária > P95. Para cada dia do
evento, o gráfico usa o menor Mmax entre as observações diurnas (06--17 h),
isto é, a condição mais restritiva do dia no cenário escolhido.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from datetime import date, datetime, timedelta
from math import ceil, floor, isfinite
from pathlib import Path


BASE = Path(__file__).resolve().parent
MODELO = BASE / "outputs"
RESULTADOS = BASE / "resultados" / "mmax_eventos"
RESULTADOS_ADULTO = RESULTADOS / "adulto_jovem"
RESULTADOS_IDOSOS = RESULTADOS / "idosos"
TEMPERATURA = BASE / "series_temporais_csv" / "temperatura_ar_1995_2024.csv"
UMIDADE = BASE / "series_temporais_csv" / "umidade_relativa_1995_2024.csv"
GRADE_MODELO = MODELO / "liveability_faixa_historica.csv"
ARQUIVO_DADOS = RESULTADOS_ADULTO / "iag-mmax-minimo-diario-por-evento.csv"
ARQUIVO_GRAFICO = RESULTADOS_ADULTO / "iag-mmax-minimo-diario-por-evento.svg"
ARQUIVO_FREQUENCIA = RESULTADOS_ADULTO / "iag-frequencia-atividades-diaria.csv"
ARQUIVO_GRAFICO_FREQUENCIA = RESULTADOS_ADULTO / "iag-frequencia-atividades-diaria.svg"
ARQUIVO_TENDENCIA = RESULTADOS_ADULTO / "iag-tendencia-mmax-diaria.csv"
ARQUIVO_GRAFICO_TENDENCIA = RESULTADOS_ADULTO / "iag-tendencia-mmax-diaria.svg"
ARQUIVO_DADOS_HORARIOS = RESULTADOS_ADULTO / "iag-mmax-horario-por-evento.csv"
ARQUIVO_FREQUENCIA_HORARIA = RESULTADOS_ADULTO / "iag-frequencia-atividades-horaria.csv"
ARQUIVO_GRAFICO_FREQUENCIA_HORARIA = RESULTADOS_ADULTO / "iag-frequencia-atividades-horaria.svg"
ARQUIVO_TENDENCIA_HORARIA = RESULTADOS_ADULTO / "iag-tendencia-mmax-horaria.csv"
ARQUIVO_GRAFICO_TENDENCIA_HORARIA = RESULTADOS_ADULTO / "iag-tendencia-mmax-horaria.svg"
ARQUIVO_DADOS_HORARIOS_IDOSOS = RESULTADOS_IDOSOS / "iag-mmax-horario-por-evento.csv"
ARQUIVO_FREQUENCIA_HORARIA_IDOSOS = RESULTADOS_IDOSOS / "iag-frequencia-atividades-horaria.csv"
ARQUIVO_GRAFICO_FREQUENCIA_HORARIA_IDOSOS = RESULTADOS_IDOSOS / "iag-frequencia-atividades-horaria.svg"
ARQUIVO_TENDENCIA_HORARIA_IDOSOS = RESULTADOS_IDOSOS / "iag-tendencia-mmax-horaria.csv"
ARQUIVO_GRAFICO_TENDENCIA_HORARIA_IDOSOS = RESULTADOS_IDOSOS / "iag-tendencia-mmax-horaria.svg"

CENARIO = ("Day-Outdoors", "Young_adult", "3")
CENARIO_IDOSOS = ("Day-Outdoors", "65_over", "3")
HORA_INICIO, HORA_FIM = 6, 18
DURACAO_MINIMA_DIAS = 3

# A categoria é a atividade mais intensa cujo equivalente em MET não excede
# Mmax. Os exemplos são os mesmos do tutorial original do projeto.
CATEGORIAS_ATIVIDADE = [
    (1.5, "#4c78a8", "Repouso ou estudo", "1,5 MET"),
    (2.3, "#72b7b2", "Caminhada lenta", "2,3 MET"),
    (3.3, "#54a24b", "Limpeza doméstica", "3,3 MET"),
    (4.0, "#eeca3b", "Subir escadas devagar", "4,0 MET"),
    (5.0, "#f58518", "Dança", "5,0 MET"),
    (6.0, "#e45756", "Corrida leve (≈ 6,4 km/h)", "6,0 MET"),
    (7.5, "#b279a2", "Carregar carga pesada", "7,5 MET"),
]


def ler_serie(caminho: Path) -> dict[datetime, float]:
    valores = {}
    with caminho.open(newline="", encoding="utf-8") as arquivo:
        for linha in csv.DictReader(arquivo):
            try:
                valor = float(linha["valor"])
            except (TypeError, ValueError):
                continue
            if valor != -99:
                valores[datetime.fromisoformat(linha["data_hora"])] = valor
    return valores


def percentil(valores: list[float], q: float) -> float:
    """Percentil com interpolação linear, equivalente ao padrão do pandas."""
    ordenados = sorted(valores)
    posicao = (len(ordenados) - 1) * q
    inferior, superior = floor(posicao), ceil(posicao)
    return ordenados[inferior] + (ordenados[superior] - ordenados[inferior]) * (
        posicao - inferior
    )


def arredondar_meio(valor: float) -> float:
    # As umidades são não negativas; evita a regra de arredondamento bancário.
    return floor(valor * 2 + 0.5) / 2


def ler_lookup_mmax(
    caminho: Path, cenario: tuple[str, str, str] = CENARIO
) -> dict[tuple[float, float], float]:
    lookup = {}
    with caminho.open(newline="", encoding="utf-8") as arquivo:
        for linha in csv.DictReader(arquivo):
            if (
                linha["exposure"], linha["profile_id"], linha["duration_h"]
            ) != cenario:
                continue
            try:
                chave = (
                    round(float(linha["air_temperature_C"]), 1),
                    round(float(linha["relative_humidity_percent"]), 1),
                )
                lookup[chave] = float(linha["Mmax_MET"])
            except (TypeError, ValueError):
                continue
    return lookup


def eventos_quentes(
    maximas_diarias: dict[date, float], limiar: float
) -> list[list[date]]:
    eventos, atual = [], []
    anterior = None
    for dia in sorted(maximas_diarias):
        acima = maximas_diarias[dia] > limiar
        if acima:
            if anterior is not None and dia == anterior + timedelta(days=1):
                atual.append(dia)
            else:
                if atual:
                    eventos.append(atual)
                atual = [dia]
            anterior = dia
        else:
            if atual:
                eventos.append(atual)
                atual = []
            anterior = None
    if atual:
        eventos.append(atual)
    return [evento for evento in eventos if len(evento) >= DURACAO_MINIMA_DIAS]


def cor_evento(indice: int, total: int) -> str:
    # Matiz distribuída em todo o círculo; cada evento recebe uma cor estável.
    return f"hsl({(indice * 137.508) % 360:.1f}, 67%, 40%)"


def escapar(texto: str) -> str:
    return texto.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def categoria_atividade(mmax: float) -> tuple[float, str, str, str]:
    categoria = CATEGORIAS_ATIVIDADE[0]
    for candidata in CATEGORIAS_ATIVIDADE:
        if mmax >= candidata[0]:
            categoria = candidata
        else:
            break
    return categoria


def gerar_svg(
    registros: list[dict], p95: float, inicio: date, fim: date,
    estacao: str = "Estação IAG", periodo: str = "1995–2024", perfil: str = "adulto jovem"
) -> str:
    largura, altura = 1920, 1010
    esquerda, direita, topo, base = 115, 420, 95, 245
    pw, ph = largura - esquerda - direita, altura - topo - base
    valores = [r["mmax_min_met"] for r in registros if isfinite(r["mmax_min_met"])]
    # Mantém visíveis os níveis de atividade de referência, inclusive os que
    # ficam acima dos valores observados nos eventos.
    minimo = min(1.5, floor((min(valores) - 0.15) * 2) / 2)
    maximo = max(7.5, ceil((max(valores) + 0.15) * 2) / 2)

    dias_totais = (fim - inicio).days
    def x(dia: date) -> float:
        return esquerda + ((dia - inicio).days / dias_totais) * pw
    def y(valor: float) -> float:
        return topo + (maximo - valor) / (maximo - minimo) * ph

    partes = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{largura}" height="{altura}" viewBox="0 0 {largura} {altura}" role="img" aria-labelledby="titulo descricao">',
        "<style>"
        "text{font-family:Arial,Helvetica,sans-serif;fill:#202124}"
        ".titulo{font-size:27px;font-weight:700}.subtitulo{font-size:16px;fill:#4b5563}"
        ".eixo{font-size:14px;fill:#374151}.legenda{font-size:12px;fill:#374151}.atividade{font-size:14px;fill:#202124}.atividade-menor{font-size:12px;fill:#4b5563}"
        ".grade{stroke:#d1d5db;stroke-width:1}.borda{stroke:#4b5563;stroke-width:1.2}"
        ".linha{fill:none;stroke-width:2.8;stroke-linejoin:round;stroke-linecap:round}"
        ".ponto{stroke:#fff;stroke-width:0.8}"
        "</style>",
        "<title id=\"titulo\">Mmax nos eventos de Tmax diária acima do percentil 95</title>",
        f"<desc id=\"descricao\">Série diária do menor Mmax diurno por evento. P95 da Tmax diária: {p95:.1f} °C.</desc>",
        '<rect width="100%" height="100%" fill="white"/>',
        '<text class="titulo" x="115" y="38">Mmax durante eventos de Tmax diária &gt; P95</text>',
        f'<text class="subtitulo" x="115" y="64">{escapar(estacao)}, {escapar(periodo)} · P95 = {p95:.1f} °C · Day-Outdoors / {escapar(perfil)} / 3 h · mínimo diário entre 06–17 h</text>',
    ]

    atividades = CATEGORIAS_ATIVIDADE

    passo = 0.5
    marca = minimo
    while marca <= maximo + 1e-9:
        yy = y(marca)
        partes.append(f'<line class="grade" x1="{esquerda}" y1="{yy:.1f}" x2="{largura-direita}" y2="{yy:.1f}"/>')
        partes.append(f'<text class="eixo" text-anchor="end" x="{esquerda-12}" y="{yy+5:.1f}">{marca:.1f}</text>')
        marca += passo
    partes.append(f'<line class="borda" x1="{esquerda}" y1="{topo}" x2="{esquerda}" y2="{topo+ph}"/>')
    partes.append(f'<line class="borda" x1="{esquerda}" y1="{topo+ph}" x2="{largura-direita}" y2="{topo+ph}"/>')
    partes.append(f'<text class="eixo" transform="translate(27 {topo+ph/2:.1f}) rotate(-90)" text-anchor="middle">Mmax (MET)</text>')

    # Linhas discretas conectam a série às atividades de referência à direita.
    for met, cor, _, _ in atividades:
        yy = y(met)
        partes.append(
            f'<line x1="{esquerda}" y1="{yy:.1f}" x2="{largura-direita}" y2="{yy:.1f}" stroke="{cor}" stroke-width="1.3" stroke-dasharray="5 5" opacity="0.30"/>'
        )

    for ano in range(inicio.year, fim.year + 1, 2):
        xx = x(date(ano, 1, 1))
        partes.append(f'<line class="grade" x1="{xx:.1f}" y1="{topo}" x2="{xx:.1f}" y2="{topo+ph}"/>')
        partes.append(f'<text class="eixo" text-anchor="middle" x="{xx:.1f}" y="{topo+ph+24}">{ano}</text>')

    por_evento = defaultdict(list)
    for registro in registros:
        por_evento[registro["evento"]].append(registro)
    for evento, pontos in por_evento.items():
        cor = cor_evento(evento - 1, len(por_evento))
        caminho = " ".join(
            f"{'M' if i == 0 else 'L'} {x(p['data']):.1f} {y(p['mmax_min_met']):.1f}"
            for i, p in enumerate(pontos)
        )
        partes.append(f'<path class="linha" d="{caminho}" stroke="{cor}"/>')
        for ponto in pontos:
            partes.append(
                f'<circle class="ponto" cx="{x(ponto["data"]):.1f}" cy="{y(ponto["mmax_min_met"]):.1f}" r="3.5" fill="{cor}"/>'
            )

    partes.append(f'<text class="subtitulo" x="{esquerda}" y="{topo+ph+55}">Cada cor identifica um período contínuo (≥ 3 dias); não há linha entre períodos diferentes.</text>')

    painel_x = largura - direita + 38
    partes.append(f'<text class="titulo" x="{painel_x}" y="38">Atividade sustentável</text>')
    partes.append(f'<text class="atividade-menor" x="{painel_x}" y="62">quando Mmax ≥ o valor indicado</text>')
    for met, cor, atividade, rotulo_met in atividades:
        yy = y(met)
        partes.append(f'<line x1="{painel_x}" y1="{yy:.1f}" x2="{painel_x+22}" y2="{yy:.1f}" stroke="{cor}" stroke-width="5"/>')
        partes.append(f'<text class="atividade" x="{painel_x+32}" y="{yy-3:.1f}">{escapar(atividade)}</text>')
        partes.append(f'<text class="atividade-menor" x="{painel_x+32}" y="{yy+14:.1f}">{rotulo_met}</text>')
    partes.append(f'<text class="atividade-menor" x="{painel_x}" y="{topo+ph+50}">Exemplos de intensidade metabólica;</text>')
    partes.append(f'<text class="atividade-menor" x="{painel_x}" y="{topo+ph+66}">não substituem avaliação clínica.</text>')

    y_legenda = topo + ph + 84
    for i, evento in enumerate(sorted(por_evento)):
        referencia = por_evento[evento][0]
        coluna, linha = divmod(i, 12)
        xx = esquerda + coluna * 420
        yy = y_legenda + linha * 20
        rotulo = f"E{evento:02d}  {referencia['inicio'].strftime('%d/%m/%Y')}–{referencia['fim'].strftime('%d/%m/%Y')}"
        partes.append(f'<line x1="{xx}" y1="{yy-4}" x2="{xx+18}" y2="{yy-4}" stroke="{cor_evento(evento-1, len(por_evento))}" stroke-width="3"/>')
        partes.append(f'<text class="legenda" x="{xx+24}" y="{yy}">{escapar(rotulo)}</text>')
    partes.append("</svg>")
    return "\n".join(partes)


def gerar_svg_frequencia(
    frequencias: list[dict], unidade: str = "dias", perfil: str = "adulto jovem",
    numero_eventos: int = 48, periodo: str = "1995–2024", estacao: str = "Estação IAG"
) -> str:
    largura, altura = 1100, 590
    esquerda, direita, topo, base = 330, 75, 110, 75
    pw, ph = largura - esquerda - direita, altura - topo - base
    maximo = max(item["dias"] for item in frequencias)
    passo_eixo = 20 if maximo <= 160 else 200
    maximo_eixo = int(ceil(maximo / passo_eixo) * passo_eixo)
    altura_barra = 32
    espacamento = ph / len(frequencias)

    partes = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{largura}" height="{altura}" viewBox="0 0 {largura} {altura}" role="img" aria-labelledby="titulo descricao">',
        "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}.titulo{font-size:27px;font-weight:700}.subtitulo{font-size:16px;fill:#4b5563}.eixo{font-size:14px;fill:#374151}.categoria{font-size:16px}.valor{font-size:16px;font-weight:700}</style>",
        f'<title id="titulo">Frequência de {unidade} por categoria de atividade sustentável</title>',
        f'<desc id="descricao">Contagem das {unidade} em eventos de Tmax acima do percentil 95, classificadas pela atividade mais intensa sustentável conforme Mmax.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text class="titulo" x="330" y="42">Frequência por atividade máxima sustentável</text>',
        f'<text class="subtitulo" x="330" y="68">{escapar(estacao)} · {unidade.capitalize()} dos {numero_eventos} eventos de Tmax diária &gt; P95 · {escapar(perfil)} ({escapar(periodo)})</text>',
    ]

    for marca in range(0, maximo_eixo + 1, passo_eixo):
        xx = esquerda + pw * marca / maximo_eixo
        partes.append(f'<line x1="{xx:.1f}" y1="{topo}" x2="{xx:.1f}" y2="{topo+ph}" stroke="#d1d5db" stroke-width="1"/>')
        partes.append(f'<text class="eixo" text-anchor="middle" x="{xx:.1f}" y="{topo+ph+25}">{marca}</text>')
    partes.append(f'<line x1="{esquerda}" y1="{topo+ph}" x2="{largura-direita}" y2="{topo+ph}" stroke="#4b5563" stroke-width="1.2"/>')
    partes.append(f'<text class="eixo" text-anchor="middle" x="{esquerda+pw/2:.1f}" y="{altura-18}">Número de {unidade}</text>')

    for indice, item in enumerate(frequencias):
        yy = topo + indice * espacamento + (espacamento - altura_barra) / 2
        largura_barra = pw * item["dias"] / maximo_eixo
        meio = yy + altura_barra / 2
        partes.append(f'<text class="categoria" text-anchor="end" x="{esquerda-18}" y="{meio-2:.1f}">{escapar(item["atividade"])}</text>')
        partes.append(f'<text class="eixo" text-anchor="end" x="{esquerda-18}" y="{meio+16:.1f}">{escapar(item["faixa"])}</text>')
        partes.append(f'<rect x="{esquerda}" y="{yy:.1f}" width="{largura_barra:.1f}" height="{altura_barra}" fill="{item["cor"]}"/>')
        partes.append(f'<text class="valor" x="{esquerda+largura_barra+10:.1f}" y="{meio+6:.1f}">{item["dias"]}</text>')
    partes.append(f'<text class="subtitulo" x="330" y="94">Cada {unidade[:-1] if unidade.endswith("s") else unidade} entra apenas na categoria mais intensa cujo limiar em MET é ≤ Mmax.</text>')
    partes.append("</svg>")
    return "\n".join(partes)


def tendencia_linear(dados: list[dict]) -> tuple[float, float, float]:
    """Retorna inclinação anual, intercepto e R² da regressão das médias anuais."""
    xs = [item["ano"] for item in dados]
    ys = [item["mmax_medio_met"] for item in dados]
    x_medio = sum(xs) / len(xs)
    y_medio = sum(ys) / len(ys)
    variacao_x = sum((x - x_medio) ** 2 for x in xs)
    inclinacao = sum((x - x_medio) * (y - y_medio) for x, y in zip(xs, ys)) / variacao_x
    intercepto = y_medio - inclinacao * x_medio
    estimados = [intercepto + inclinacao * x for x in xs]
    soma_residuos = sum((y - previsto) ** 2 for y, previsto in zip(ys, estimados))
    soma_total = sum((y - y_medio) ** 2 for y in ys)
    r2 = 1 - soma_residuos / soma_total if soma_total else 0.0
    return inclinacao, intercepto, r2


def gerar_svg_tendencia(
    dados: list[dict], inclinacao: float, intercepto: float, r2: float,
    unidade: str = "dias", perfil: str = "adulto jovem", ano_inicio: int = 1995,
    ano_fim: int = 2024, estacao: str = "Estação IAG"
) -> str:
    largura, altura = 1200, 690
    esquerda, direita, topo, base = 115, 55, 115, 120
    pw, ph = largura - esquerda - direita, altura - topo - base
    anos = [item["ano"] for item in dados]
    valores = [item["mmax_medio_met"] for item in dados]
    minimo = floor((min(valores) - 0.15) * 2) / 2
    maximo = ceil((max(valores) + 0.15) * 2) / 2
    if minimo == maximo:
        maximo += 0.5

    def x(ano: int) -> float:
        return esquerda + (ano - ano_inicio) / (ano_fim - ano_inicio) * pw
    def y(valor: float) -> float:
        return topo + (maximo - valor) / (maximo - minimo) * ph

    sinal = "+" if inclinacao >= 0 else "−"
    partes = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{largura}" height="{altura}" viewBox="0 0 {largura} {altura}" role="img" aria-labelledby="titulo descricao">',
        "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}.titulo{font-size:27px;font-weight:700}.subtitulo{font-size:16px;fill:#4b5563}.eixo{font-size:14px;fill:#374151}.valor{font-size:13px;fill:#374151}.nota{font-size:14px;fill:#4b5563}</style>",
        f'<title id="titulo">Tendência anual de Mmax em eventos quentes, baseada em {unidade}</title>',
        f'<desc id="descricao">Média anual de Mmax em {unidade} de eventos de Tmax acima do percentil 95, com linha de regressão linear.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text class="titulo" x="115" y="42">Tendência anual de Mmax em eventos quentes</text>',
        f'<text class="subtitulo" x="115" y="68">{escapar(estacao)} · Média anual de Mmax em {unidade} diurnas (06–17 h) · Day-Outdoors / {escapar(perfil)} / 3 h</text>',
        f'<text class="subtitulo" x="115" y="92">Tendência linear: {sinal}{abs(inclinacao * 10):.2f} MET por década · R² = {r2:.2f} · {len(dados)} anos com eventos</text>',
    ]
    marca = minimo
    while marca <= maximo + 1e-9:
        yy = y(marca)
        partes.append(f'<line x1="{esquerda}" y1="{yy:.1f}" x2="{largura-direita}" y2="{yy:.1f}" stroke="#d1d5db" stroke-width="1"/>')
        partes.append(f'<text class="eixo" text-anchor="end" x="{esquerda-12}" y="{yy+5:.1f}">{marca:.1f}</text>')
        marca += 0.5
    for ano in range(ano_inicio, ano_fim + 1, 5):
        xx = x(ano)
        partes.append(f'<line x1="{xx:.1f}" y1="{topo}" x2="{xx:.1f}" y2="{topo+ph}" stroke="#e5e7eb" stroke-width="1"/>')
        partes.append(f'<text class="eixo" text-anchor="middle" x="{xx:.1f}" y="{topo+ph+25}">{ano}</text>')
    partes.append(f'<line x1="{esquerda}" y1="{topo+ph}" x2="{largura-direita}" y2="{topo+ph}" stroke="#4b5563" stroke-width="1.2"/>')
    partes.append(f'<line x1="{esquerda}" y1="{topo}" x2="{esquerda}" y2="{topo+ph}" stroke="#4b5563" stroke-width="1.2"/>')
    partes.append(f'<text class="eixo" transform="translate(27 {topo+ph/2:.1f}) rotate(-90)" text-anchor="middle">Mmax médio (MET)</text>')
    partes.append(f'<text class="eixo" text-anchor="middle" x="{esquerda+pw/2:.1f}" y="{altura-25}">Ano</text>')

    linha_inicio = intercepto + inclinacao * ano_inicio
    linha_fim = intercepto + inclinacao * ano_fim
    partes.append(f'<path d="M {x(ano_inicio):.1f} {y(linha_inicio):.1f} L {x(ano_fim):.1f} {y(linha_fim):.1f}" fill="none" stroke="#c43c39" stroke-width="3"/>')
    for item in dados:
        xx, yy = x(item["ano"]), y(item["mmax_medio_met"])
        partes.append(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="5" fill="#2f5597" stroke="white" stroke-width="1.2"/>')
        partes.append(f'<text class="valor" text-anchor="middle" x="{xx:.1f}" y="{yy-11:.1f}">n={item["dias"]}</text>')
    partes.append('<line x1="115" y1="620" x2="140" y2="620" stroke="#2f5597" stroke-width="3"/><text class="nota" x="150" y="625">Média anual</text>')
    partes.append('<line x1="310" y1="620" x2="335" y2="620" stroke="#c43c39" stroke-width="3"/><text class="nota" x="345" y="625">Regressão linear</text>')
    partes.append(f'<text class="nota" x="115" y="655">Análise condicionada às {unidade} de eventos de Tmax &gt; P95; anos sem eventos não entram no ajuste.</text>')
    partes.append("</svg>")
    return "\n".join(partes)


def main() -> None:
    temperatura = ler_serie(TEMPERATURA)
    umidade = ler_serie(UMIDADE)
    por_dia = defaultdict(list)
    for instante, valor in temperatura.items():
        por_dia[instante.date()].append(valor)
    maximas_diarias = {dia: max(valores) for dia, valores in por_dia.items()}
    p95 = percentil(list(maximas_diarias.values()), 0.95)
    eventos = eventos_quentes(maximas_diarias, p95)
    lookup = ler_lookup_mmax(GRADE_MODELO)

    registros = []
    for numero, evento in enumerate(eventos, start=1):
        for dia in evento:
            mmax_do_dia = []
            for hora in range(HORA_INICIO, HORA_FIM):
                instante = datetime.combine(dia, datetime.min.time()) + timedelta(hours=hora)
                if instante not in temperatura or instante not in umidade:
                    continue
                chave = (round(temperatura[instante], 1), arredondar_meio(umidade[instante]))
                mmax = lookup.get(chave)
                if mmax is not None and isfinite(mmax):
                    mmax_do_dia.append(mmax)
            if mmax_do_dia:
                registros.append({
                    "evento": numero,
                    "inicio": evento[0],
                    "fim": evento[-1],
                    "data": dia,
                    "duracao_dias": len(evento),
                    "tmax_diaria_c": maximas_diarias[dia],
                    "mmax_min_met": min(mmax_do_dia),
                    "observacoes_diurnas": len(mmax_do_dia),
                })

    RESULTADOS.mkdir(parents=True, exist_ok=True)
    RESULTADOS_ADULTO.mkdir(exist_ok=True)
    RESULTADOS_IDOSOS.mkdir(exist_ok=True)
    with ARQUIVO_DADOS.open("w", newline="", encoding="utf-8") as arquivo:
        campos = list(registros[0])
        escrita = csv.DictWriter(arquivo, fieldnames=campos)
        escrita.writeheader()
        for registro in registros:
            escrita.writerow({chave: valor.isoformat() if isinstance(valor, date) else valor for chave, valor in registro.items()})

    inicio, fim = min(maximas_diarias), max(maximas_diarias)
    ARQUIVO_GRAFICO.write_text(gerar_svg(registros, p95, inicio, fim), encoding="utf-8")

    contagens = {categoria[2]: 0 for categoria in CATEGORIAS_ATIVIDADE}
    for registro in registros:
        contagens[categoria_atividade(registro["mmax_min_met"])[2]] += 1
    frequencias = []
    for indice, (met, cor, atividade, rotulo_met) in enumerate(CATEGORIAS_ATIVIDADE):
        proximo = CATEGORIAS_ATIVIDADE[indice + 1][0] if indice + 1 < len(CATEGORIAS_ATIVIDADE) else None
        faixa = f"{rotulo_met} a < {proximo:.1f} MET" if proximo else f"≥ {rotulo_met}"
        frequencias.append({"atividade": atividade, "faixa": faixa.replace(".", ","), "dias": contagens[atividade], "cor": cor})
    with ARQUIVO_FREQUENCIA.open("w", newline="", encoding="utf-8") as arquivo:
        escrita = csv.DictWriter(arquivo, fieldnames=["atividade", "faixa", "dias"])
        escrita.writeheader()
        escrita.writerows({chave: item[chave] for chave in ["atividade", "faixa", "dias"]} for item in frequencias)
    ARQUIVO_GRAFICO_FREQUENCIA.write_text(gerar_svg_frequencia(frequencias), encoding="utf-8")

    por_ano = defaultdict(list)
    for registro in registros:
        por_ano[registro["data"].year].append(registro["mmax_min_met"])
    tendencia_anual = [
        {
            "ano": ano,
            "dias": len(valores),
            "mmax_medio_met": sum(valores) / len(valores),
            "mmax_minimo_met": min(valores),
            "mmax_maximo_met": max(valores),
        }
        for ano, valores in sorted(por_ano.items())
    ]
    inclinacao, intercepto, r2 = tendencia_linear(tendencia_anual)
    with ARQUIVO_TENDENCIA.open("w", newline="", encoding="utf-8") as arquivo:
        escrita = csv.DictWriter(arquivo, fieldnames=list(tendencia_anual[0]))
        escrita.writeheader()
        escrita.writerows(tendencia_anual)
    ARQUIVO_GRAFICO_TENDENCIA.write_text(
        gerar_svg_tendencia(tendencia_anual, inclinacao, intercepto, r2), encoding="utf-8"
    )

    # Análise horária: preserva cada Mmax do período quente, sem reduzi-lo ao
    # pior valor de cada dia. Como o cenário é Day-Outdoors, entram 06--17 h.
    evento_por_dia = {
        dia: (numero, evento[0], evento[-1], len(evento))
        for numero, evento in enumerate(eventos, start=1)
        for dia in evento
    }
    registros_horarios = []
    for instante in sorted(temperatura):
        if instante.date() not in evento_por_dia or not HORA_INICIO <= instante.hour < HORA_FIM:
            continue
        if instante not in umidade:
            continue
        chave = (round(temperatura[instante], 1), arredondar_meio(umidade[instante]))
        mmax = lookup.get(chave)
        if mmax is None or not isfinite(mmax):
            continue
        numero, inicio_evento, fim_evento, duracao = evento_por_dia[instante.date()]
        registros_horarios.append({
            "evento": numero,
            "inicio": inicio_evento,
            "fim": fim_evento,
            "duracao_dias": duracao,
            "data_hora": instante,
            "temperatura_ar_c": temperatura[instante],
            "umidade_relativa_percent": umidade[instante],
            "mmax_met": mmax,
        })
    with ARQUIVO_DADOS_HORARIOS.open("w", newline="", encoding="utf-8") as arquivo:
        campos = list(registros_horarios[0])
        escrita = csv.DictWriter(arquivo, fieldnames=campos)
        escrita.writeheader()
        for registro in registros_horarios:
            escrita.writerow({
                chave: valor.isoformat(sep=" ") if isinstance(valor, datetime)
                else valor.isoformat() if isinstance(valor, date)
                else valor
                for chave, valor in registro.items()
            })

    contagens_horarias = {categoria[2]: 0 for categoria in CATEGORIAS_ATIVIDADE}
    for registro in registros_horarios:
        contagens_horarias[categoria_atividade(registro["mmax_met"])[2]] += 1
    frequencias_horarias = []
    for indice, (met, cor, atividade, rotulo_met) in enumerate(CATEGORIAS_ATIVIDADE):
        proximo = CATEGORIAS_ATIVIDADE[indice + 1][0] if indice + 1 < len(CATEGORIAS_ATIVIDADE) else None
        faixa = f"{rotulo_met} a < {proximo:.1f} MET" if proximo else f"≥ {rotulo_met}"
        frequencias_horarias.append({
            "atividade": atividade,
            "faixa": faixa.replace(".", ","),
            "dias": contagens_horarias[atividade],
            "cor": cor,
        })
    with ARQUIVO_FREQUENCIA_HORARIA.open("w", newline="", encoding="utf-8") as arquivo:
        escrita = csv.DictWriter(arquivo, fieldnames=["atividade", "faixa", "horas"])
        escrita.writeheader()
        for item in frequencias_horarias:
            escrita.writerow({"atividade": item["atividade"], "faixa": item["faixa"], "horas": item["dias"]})
    ARQUIVO_GRAFICO_FREQUENCIA_HORARIA.write_text(
        gerar_svg_frequencia(frequencias_horarias, unidade="horas"), encoding="utf-8"
    )

    por_ano_horario = defaultdict(list)
    for registro in registros_horarios:
        por_ano_horario[registro["data_hora"].year].append(registro["mmax_met"])
    tendencia_horaria = [
        {
            "ano": ano,
            "dias": len(valores),  # chave interna reutilizada pelo gerador de SVG
            "mmax_medio_met": sum(valores) / len(valores),
            "mmax_minimo_met": min(valores),
            "mmax_maximo_met": max(valores),
        }
        for ano, valores in sorted(por_ano_horario.items())
    ]
    inclinacao_horaria, intercepto_horario, r2_horario = tendencia_linear(tendencia_horaria)
    with ARQUIVO_TENDENCIA_HORARIA.open("w", newline="", encoding="utf-8") as arquivo:
        escrita = csv.DictWriter(
            arquivo,
            fieldnames=["ano", "horas", "mmax_medio_met", "mmax_minimo_met", "mmax_maximo_met"],
        )
        escrita.writeheader()
        for item in tendencia_horaria:
            escrita.writerow({
                "ano": item["ano"], "horas": item["dias"],
                "mmax_medio_met": item["mmax_medio_met"],
                "mmax_minimo_met": item["mmax_minimo_met"],
                "mmax_maximo_met": item["mmax_maximo_met"],
            })
    ARQUIVO_GRAFICO_TENDENCIA_HORARIA.write_text(
        gerar_svg_tendencia(
            tendencia_horaria, inclinacao_horaria, intercepto_horario, r2_horario, unidade="horas"
        ),
        encoding="utf-8",
    )

    # Repete a análise horária para o perfil 65_over, mantendo os mesmos
    # períodos de calor e horários de exposição do cenário externo.
    lookup_idosos = ler_lookup_mmax(GRADE_MODELO, CENARIO_IDOSOS)
    registros_horarios_idosos = []
    for instante in sorted(temperatura):
        if instante.date() not in evento_por_dia or not HORA_INICIO <= instante.hour < HORA_FIM:
            continue
        if instante not in umidade:
            continue
        chave = (round(temperatura[instante], 1), arredondar_meio(umidade[instante]))
        mmax = lookup_idosos.get(chave)
        if mmax is None or not isfinite(mmax):
            continue
        numero, inicio_evento, fim_evento, duracao = evento_por_dia[instante.date()]
        registros_horarios_idosos.append({
            "evento": numero,
            "inicio": inicio_evento,
            "fim": fim_evento,
            "duracao_dias": duracao,
            "data_hora": instante,
            "temperatura_ar_c": temperatura[instante],
            "umidade_relativa_percent": umidade[instante],
            "mmax_met": mmax,
        })
    with ARQUIVO_DADOS_HORARIOS_IDOSOS.open("w", newline="", encoding="utf-8") as arquivo:
        campos = list(registros_horarios_idosos[0])
        escrita = csv.DictWriter(arquivo, fieldnames=campos)
        escrita.writeheader()
        for registro in registros_horarios_idosos:
            escrita.writerow({
                chave: valor.isoformat(sep=" ") if isinstance(valor, datetime)
                else valor.isoformat() if isinstance(valor, date)
                else valor
                for chave, valor in registro.items()
            })

    contagens_idosos = {categoria[2]: 0 for categoria in CATEGORIAS_ATIVIDADE}
    for registro in registros_horarios_idosos:
        contagens_idosos[categoria_atividade(registro["mmax_met"])[2]] += 1
    frequencias_idosos = []
    for indice, (met, cor, atividade, rotulo_met) in enumerate(CATEGORIAS_ATIVIDADE):
        proximo = CATEGORIAS_ATIVIDADE[indice + 1][0] if indice + 1 < len(CATEGORIAS_ATIVIDADE) else None
        faixa = f"{rotulo_met} a < {proximo:.1f} MET" if proximo else f"≥ {rotulo_met}"
        frequencias_idosos.append({
            "atividade": atividade,
            "faixa": faixa.replace(".", ","),
            "dias": contagens_idosos[atividade],
            "cor": cor,
        })
    with ARQUIVO_FREQUENCIA_HORARIA_IDOSOS.open("w", newline="", encoding="utf-8") as arquivo:
        escrita = csv.DictWriter(arquivo, fieldnames=["atividade", "faixa", "horas"])
        escrita.writeheader()
        for item in frequencias_idosos:
            escrita.writerow({"atividade": item["atividade"], "faixa": item["faixa"], "horas": item["dias"]})
    ARQUIVO_GRAFICO_FREQUENCIA_HORARIA_IDOSOS.write_text(
        gerar_svg_frequencia(frequencias_idosos, unidade="horas", perfil="65 anos ou mais"),
        encoding="utf-8",
    )

    por_ano_idosos = defaultdict(list)
    for registro in registros_horarios_idosos:
        por_ano_idosos[registro["data_hora"].year].append(registro["mmax_met"])
    tendencia_idosos = [
        {
            "ano": ano,
            "dias": len(valores),
            "mmax_medio_met": sum(valores) / len(valores),
            "mmax_minimo_met": min(valores),
            "mmax_maximo_met": max(valores),
        }
        for ano, valores in sorted(por_ano_idosos.items())
    ]
    inclinacao_idosos, intercepto_idosos, r2_idosos = tendencia_linear(tendencia_idosos)
    with ARQUIVO_TENDENCIA_HORARIA_IDOSOS.open("w", newline="", encoding="utf-8") as arquivo:
        escrita = csv.DictWriter(
            arquivo,
            fieldnames=["ano", "horas", "mmax_medio_met", "mmax_minimo_met", "mmax_maximo_met"],
        )
        escrita.writeheader()
        for item in tendencia_idosos:
            escrita.writerow({
                "ano": item["ano"], "horas": item["dias"],
                "mmax_medio_met": item["mmax_medio_met"],
                "mmax_minimo_met": item["mmax_minimo_met"],
                "mmax_maximo_met": item["mmax_maximo_met"],
            })
    ARQUIVO_GRAFICO_TENDENCIA_HORARIA_IDOSOS.write_text(
        gerar_svg_tendencia(
            tendencia_idosos, inclinacao_idosos, intercepto_idosos, r2_idosos,
            unidade="horas", perfil="65 anos ou mais"
        ),
        encoding="utf-8",
    )

    print(f"P95 da Tmax diária: {p95:.1f} °C")
    print(f"Eventos (≥ {DURACAO_MINIMA_DIAS} dias): {len(eventos)}")
    print(f"Dados: {ARQUIVO_DADOS}")
    print(f"Gráfico: {ARQUIVO_GRAFICO}")
    print(f"Frequências: {ARQUIVO_FREQUENCIA}")
    print(f"Gráfico de frequências: {ARQUIVO_GRAFICO_FREQUENCIA}")
    print(f"Tendência: {ARQUIVO_TENDENCIA}")
    print(f"Gráfico de tendência: {ARQUIVO_GRAFICO_TENDENCIA}")
    print(f"Inclinação: {inclinacao * 10:+.3f} MET por década; R² = {r2:.3f}")
    print(f"Dados horários: {ARQUIVO_DADOS_HORARIOS}")
    print(f"Frequência horária: {ARQUIVO_FREQUENCIA_HORARIA}")
    print(f"Gráfico de frequência horária: {ARQUIVO_GRAFICO_FREQUENCIA_HORARIA}")
    print(f"Tendência horária: {ARQUIVO_TENDENCIA_HORARIA}")
    print(f"Gráfico de tendência horária: {ARQUIVO_GRAFICO_TENDENCIA_HORARIA}")
    print(f"Inclinação horária: {inclinacao_horaria * 10:+.3f} MET por década; R² = {r2_horario:.3f}")
    print(f"Dados horários (idosos): {ARQUIVO_DADOS_HORARIOS_IDOSOS}")
    print(f"Frequência horária (idosos): {ARQUIVO_FREQUENCIA_HORARIA_IDOSOS}")
    print(f"Gráfico de frequência horária (idosos): {ARQUIVO_GRAFICO_FREQUENCIA_HORARIA_IDOSOS}")
    print(f"Tendência horária (idosos): {ARQUIVO_TENDENCIA_HORARIA_IDOSOS}")
    print(f"Gráfico de tendência horária (idosos): {ARQUIVO_GRAFICO_TENDENCIA_HORARIA_IDOSOS}")
    print(f"Inclinação horária (idosos): {inclinacao_idosos * 10:+.3f} MET por década; R² = {r2_idosos:.3f}")


if __name__ == "__main__":
    main()
