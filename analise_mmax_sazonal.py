"""Mmax médio sazonal e tendências: IAG e Mirante de Santana.

Esta análise usa todas as observações diurnas locais (06--17 h) disponíveis,
avaliadas no cenário Day-Outdoors de 3 h. As estações são meteorológicas do
Hemisfério Sul: verão (dez--fev), outono (mar--mai), inverno (jun--ago) e
primavera (set--nov). Dezembro compõe o verão do ano seguinte.
"""

from __future__ import annotations

import csv
import importlib.util
from collections import defaultdict
from datetime import datetime
from math import ceil, floor, isfinite
from pathlib import Path


RAIZ = Path(__file__).resolve().parent
SCRIPT_IAG = RAIZ / "dados_estacao_iag" / "grafico_mmax_eventos_tmax_p95.py"
ARQUIVO_MODELO = RAIZ / "dados_estacao_iag" / "outputs" / "liveability_faixa_historica.csv"
SAIDA = RAIZ / "analises_comparativas" / "mmax_sazonal"
HORA_INICIO, HORA_FIM = 6, 18

PERFIS = {
    "adulto_jovem": (("Day-Outdoors", "Young_adult", "3"), "Adulto jovem"),
    "idosos": (("Day-Outdoors", "65_over", "3"), "65 anos ou mais"),
}
ESTACOES = {
    "IAG": "Estação IAG",
    "santana": "Mirante de Santana (INMET A701)",
}
ESTACOES_ANO = ("Verão", "Outono", "Inverno", "Primavera")
CORES = {"Verão": "#e45756", "Outono": "#f2cf5b", "Inverno": "#4c78a8", "Primavera": "#54a24b"}


def carregar_iag():
    spec = importlib.util.spec_from_file_location("analise_iag", SCRIPT_IAG)
    modulo = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(modulo)
    return modulo


IAG = carregar_iag()


def ler_serie_iag(caminho: Path) -> dict[datetime, float]:
    serie = {}
    with caminho.open(newline="", encoding="utf-8") as fonte:
        for linha in csv.DictReader(fonte):
            try:
                valor = float(linha["valor"])
                if valor == -99:
                    continue
                serie[datetime.fromisoformat(linha["data_hora"])] = valor
            except (KeyError, ValueError):
                continue
    return serie


def observacoes_iag() -> dict[datetime, tuple[float, float]]:
    base = RAIZ / "dados_estacao_iag" / "series_temporais_csv"
    temperatura = ler_serie_iag(base / "temperatura_ar_1995_2024.csv")
    umidade = ler_serie_iag(base / "umidade_relativa_1995_2024.csv")
    return {instante: (temperatura[instante], umidade[instante])
            for instante in temperatura.keys() & umidade.keys()}


def observacoes_santana() -> dict[datetime, tuple[float, float]]:
    caminho = RAIZ / "dados_estacao_santana" / "dados_processados" / "santana-temperatura-umidade-horaria-2006-2025.csv"
    observacoes = {}
    with caminho.open(newline="", encoding="utf-8") as fonte:
        for linha in csv.DictReader(fonte):
            try:
                observacoes[datetime.fromisoformat(linha["data_hora"])] = (
                    float(linha["temperatura_ar_c"]), float(linha["umidade_relativa_percent"])
                )
            except (KeyError, ValueError):
                continue
    return observacoes


def chave_sazonal(instante: datetime) -> tuple[str, int]:
    mes = instante.month
    if mes in (12, 1, 2):
        return "Verão", instante.year + 1 if mes == 12 else instante.year
    if mes in (3, 4, 5):
        return "Outono", instante.year
    if mes in (6, 7, 8):
        return "Inverno", instante.year
    return "Primavera", instante.year


def meses_esperados(estacao: str, ano: int) -> set[tuple[int, int]]:
    if estacao == "Verão":
        return {(ano - 1, 12), (ano, 1), (ano, 2)}
    inicio = {"Outono": 3, "Inverno": 6, "Primavera": 9}[estacao]
    return {(ano, inicio), (ano, inicio + 1), (ano, inicio + 2)}


def tendencias(registros: list[dict]) -> list[dict]:
    saida = []
    for (estacao_id, perfil_id, estacao), grupo in agrupar(registros, ("estacao_id", "perfil_id", "estacao")):
        grupo = sorted(grupo, key=lambda item: item["ano_estacao"])
        dados = [{"ano": item["ano_estacao"], "mmax_medio_met": item["mmax_medio_met"]} for item in grupo]
        inclinacao, _, r2 = IAG.tendencia_linear(dados)
        saida.append({
            "estacao_id": estacao_id,
            "estacao_meteorologica": ESTACOES[estacao_id],
            "perfil_id": perfil_id,
            "perfil": PERFIS[perfil_id][1],
            "estacao_do_ano": estacao,
            "anos_completos": len(dados),
            "tendencia_met_por_decada": inclinacao * 10,
            "r_quadrado": r2,
        })
    return saida


def agrupar(registros: list[dict], chaves: tuple[str, ...]):
    grupos = defaultdict(list)
    for registro in registros:
        grupos[tuple(registro[chave] for chave in chaves)].append(registro)
    return grupos.items()


def svg_sazonal(registros: list[dict], tendencia: list[dict], estacao_id: str, perfil_id: str) -> str:
    dados = [r for r in registros if r["estacao_id"] == estacao_id and r["perfil_id"] == perfil_id]
    linhas = [r for r in tendencia if r["estacao_id"] == estacao_id and r["perfil_id"] == perfil_id]
    largura, altura = 1440, 880
    margem_x, topo, rodape = 95, 115, 75
    painel_l, painel_t, painel_w, painel_h = 90, 110, 620, 295
    valores = [r["mmax_medio_met"] for r in dados]
    y_min = floor((min(valores) - 0.15) * 2) / 2
    y_max = ceil((max(valores) + 0.15) * 2) / 2
    anos = [r["ano_estacao"] for r in dados]
    ano_min, ano_max = min(anos), max(anos)
    if ano_min == ano_max:
        ano_max += 1

    def painel(indice: int):
        coluna, linha = indice % 2, indice // 2
        return painel_l + coluna * (painel_w + 80), painel_t + linha * (painel_h + 115)

    partes = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{largura}" height="{altura}" viewBox="0 0 {largura} {altura}" role="img" aria-labelledby="titulo descricao">',
        "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#202124}.titulo{font-size:28px;font-weight:700}.subtitulo{font-size:16px;fill:#4b5563}.painel{font-size:20px;font-weight:700}.eixo{font-size:13px;fill:#374151}.nota{font-size:13px;fill:#4b5563}</style>",
        '<title id="titulo">Mmax médio sazonal e tendências</title>',
        '<desc id="descricao">Quatro painéis com as médias anuais de Mmax para verão, outono, inverno e primavera.</desc>',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text class="titulo" x="90" y="42">Mmax médio sazonal e tendência</text>',
        f'<text class="subtitulo" x="90" y="70">{ESTACOES[estacao_id]} · {PERFIS[perfil_id][1]} · Day-Outdoors / 3 h · horas locais entre 06–17 h</text>',
        '<text class="nota" x="90" y="94">Verão: dez–fev, atribuído ao ano de jan–fev; apenas estações com os três meses disponíveis.</text>',
    ]
    for indice, estacao in enumerate(ESTACOES_ANO):
        x0, y0 = painel(indice)
        x1, y1 = x0 + painel_w, y0 + painel_h
        serie = sorted((r for r in dados if r["estacao"] == estacao), key=lambda r: r["ano_estacao"])
        linha = next(r for r in linhas if r["estacao_do_ano"] == estacao)
        def x(ano): return x0 + (ano - ano_min) / (ano_max - ano_min) * painel_w
        def y(valor): return y0 + (y_max - valor) / (y_max - y_min) * painel_h
        partes.append(f'<text class="painel" x="{x0}" y="{y0-14}">{estacao}</text>')
        marca = y_min
        while marca <= y_max + 1e-9:
            yy = y(marca)
            partes.append(f'<line x1="{x0}" y1="{yy:.1f}" x2="{x1}" y2="{yy:.1f}" stroke="#d1d5db" stroke-width="1"/>')
            partes.append(f'<text class="eixo" text-anchor="end" x="{x0-8}" y="{yy+4:.1f}">{marca:.1f}</text>')
            marca += .5
        for ano in range(((ano_min + 4) // 5) * 5, ano_max + 1, 5):
            xx = x(ano)
            partes.append(f'<line x1="{xx:.1f}" y1="{y0}" x2="{xx:.1f}" y2="{y1}" stroke="#e5e7eb" stroke-width="1"/>')
            partes.append(f'<text class="eixo" text-anchor="middle" x="{xx:.1f}" y="{y1+22}">{ano}</text>')
        partes.append(f'<line x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}" stroke="#4b5563"/>')
        partes.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" stroke="#4b5563"/>')
        inclinacao = linha["tendencia_met_por_decada"] / 10
        intercepto = sum(r["mmax_medio_met"] - inclinacao * r["ano_estacao"] for r in serie) / len(serie)
        cor = CORES[estacao]
        partes.append(f'<path d="M {x(ano_min):.1f} {y(intercepto+inclinacao*ano_min):.1f} L {x(ano_max):.1f} {y(intercepto+inclinacao*ano_max):.1f}" fill="none" stroke="#222" stroke-width="2.3"/>')
        partes.append('<path d="" fill="none"/>')
        for r in serie:
            partes.append(f'<circle cx="{x(r["ano_estacao"]):.1f}" cy="{y(r["mmax_medio_met"]):.1f}" r="4.5" fill="{cor}" stroke="white" stroke-width="1"/>')
        sinal = "+" if linha["tendencia_met_por_decada"] >= 0 else "−"
        partes.append(f'<text class="nota" x="{x0}" y="{y1+45}">Tendência: {sinal}{abs(linha["tendencia_met_por_decada"]):.2f} MET/década · R²={linha["r_quadrado"]:.2f} · n={linha["anos_completos"]} anos</text>')
    partes.append('<text class="nota" x="90" y="850">Pontos: média de todas as horas diurnas válidas em cada estação do ano. Linha preta: regressão linear das médias sazonais.</text>')
    partes.append("</svg>")
    return "\n".join(partes)


def escrever_csv(caminho: Path, registros: list[dict], campos: list[str] | None = None):
    with caminho.open("w", newline="", encoding="utf-8") as destino:
        escrita = csv.DictWriter(destino, fieldnames=campos or list(registros[0]))
        escrita.writeheader()
        escrita.writerows(registros)


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    observacoes = {"IAG": observacoes_iag(), "santana": observacoes_santana()}
    registros_sazonais = []
    cobertura = []
    for estacao_id, serie in observacoes.items():
        for perfil_id, (cenario, _) in PERFIS.items():
            lookup = IAG.ler_lookup_mmax(ARQUIVO_MODELO, cenario)
            valores, meses = defaultdict(list), defaultdict(set)
            sem_cobertura = 0
            for instante, (temperatura, umidade) in serie.items():
                if not HORA_INICIO <= instante.hour < HORA_FIM:
                    continue
                mmax = lookup.get((round(temperatura, 1), IAG.arredondar_meio(umidade)))
                if mmax is None or not isfinite(mmax):
                    sem_cobertura += 1
                    continue
                estacao, ano = chave_sazonal(instante)
                chave = (estacao, ano)
                valores[chave].append(mmax)
                meses[chave].add((instante.year, instante.month))
            for (estacao, ano), serie_sazonal in valores.items():
                if meses[(estacao, ano)] != meses_esperados(estacao, ano):
                    continue
                registros_sazonais.append({
                    "estacao_id": estacao_id,
                    "estacao_meteorologica": ESTACOES[estacao_id],
                    "perfil_id": perfil_id,
                    "perfil": PERFIS[perfil_id][1],
                    "estacao": estacao,
                    "ano_estacao": ano,
                    "horas_validas": len(serie_sazonal),
                    "mmax_medio_met": sum(serie_sazonal) / len(serie_sazonal),
                    "mmax_minimo_met": min(serie_sazonal),
                    "mmax_maximo_met": max(serie_sazonal),
                })
            cobertura.append({"estacao_id": estacao_id, "perfil_id": perfil_id, "horas_sem_cobertura_modelo": sem_cobertura})
    registros_sazonais.sort(key=lambda r: (r["estacao_id"], r["perfil_id"], r["estacao"], r["ano_estacao"]))
    resultados_tendencia = tendencias(registros_sazonais)
    escrever_csv(SAIDA / "mmax-medio-sazonal-anual-iag-e-santana.csv", registros_sazonais)
    escrever_csv(SAIDA / "tendencias-mmax-sazonal-iag-e-santana.csv", resultados_tendencia)
    escrever_csv(SAIDA / "cobertura-grade-modelo-iag-e-santana.csv", cobertura)
    for estacao_id in ESTACOES:
        for perfil_id in PERFIS:
            (SAIDA / f"mmax-medio-sazonal-{estacao_id.lower()}-{perfil_id.replace('_', '-')}.svg").write_text(
                svg_sazonal(registros_sazonais, resultados_tendencia, estacao_id, perfil_id), encoding="utf-8"
            )
    print(f"Registros sazonais completos: {len(registros_sazonais)}")
    for linha in resultados_tendencia:
        print(f"{linha['estacao_meteorologica']} | {linha['perfil']} | {linha['estacao_do_ano']}: "
              f"{linha['tendencia_met_por_decada']:+.3f} MET/década (R²={linha['r_quadrado']:.3f})")


if __name__ == "__main__":
    main()
