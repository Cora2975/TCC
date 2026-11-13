"""Análise de eventos quentes e Mmax para a estação INMET A701, Mirante.

Os arquivos brutos do INMET são horários UTC. Eles são convertidos para a hora
local de São Paulo antes do cálculo de Tmax diária, eventos e recorte diurno.
O modelo HEATLim é o mesmo usado na análise da estação IAG.
"""

from __future__ import annotations

import csv
import importlib.util
from collections import defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from zoneinfo import ZoneInfo


BASE = Path(__file__).resolve().parent
RAIZ = BASE.parent
ARQUIVO_MODELO = RAIZ / "dados_estacao_iag" / "outputs" / "liveability_faixa_historica.csv"
SCRIPT_IAG = RAIZ / "dados_estacao_iag" / "grafico_mmax_eventos_tmax_p95.py"
DADOS_PROCESSADOS = BASE / "dados_processados"
SAIDA = BASE / "resultados" / "mmax_eventos"
FUSO = ZoneInfo("America/Sao_Paulo")
HORA_INICIO, HORA_FIM = 6, 18
DURACAO_MINIMA_DIAS = 3
PERIODO = "2006–2025"

CENARIOS = {
    "adulto_jovem": (("Day-Outdoors", "Young_adult", "3"), "adulto jovem"),
    "idosos": (("Day-Outdoors", "65_over", "3"), "65 anos ou mais"),
}


def carregar_funcoes_iag():
    spec = importlib.util.spec_from_file_location("analise_iag", SCRIPT_IAG)
    modulo = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(modulo)
    return modulo


IAG = carregar_funcoes_iag()


def numero(texto: str) -> float:
    valor = float(texto.strip().replace(",", "."))
    if valor <= -9990:
        raise ValueError("valor ausente")
    return valor


def localizar_cabecalho(linhas) -> list[str]:
    for linha in linhas:
        primeira = linha[0].strip().upper() if linha else ""
        segunda = linha[1].upper() if len(linha) > 1 else ""
        if (primeira == "DATA" or primeira.startswith("DATA (")) and "HORA" in segunda:
            return linha
    raise ValueError("Cabeçalho de dados não encontrado")


def ler_observacoes() -> dict[datetime, tuple[float, float]]:
    """Retorna temperatura e UR por hora local, usando todos os CSVs A701."""
    observacoes = {}
    arquivos = sorted(BASE.glob("20[0-9][0-9]/*.CSV"))
    for arquivo in arquivos:
        with arquivo.open(encoding="latin-1", newline="") as fonte:
            linhas = csv.reader(fonte, delimiter=";")
            cabecalho = localizar_cabecalho(linhas)
            indice_temperatura = next(
                i for i, coluna in enumerate(cabecalho)
                if "TEMPERATURA DO AR - BULBO SECO" in coluna.upper()
            )
            indice_umidade = next(
                i for i, coluna in enumerate(cabecalho)
                if "UMIDADE RELATIVA DO AR, HORARIA" in coluna.upper()
            )
            for linha in linhas:
                if len(linha) <= max(indice_temperatura, indice_umidade):
                    continue
                try:
                    temperatura = numero(linha[indice_temperatura])
                    umidade = numero(linha[indice_umidade])
                    formato_data = "%Y-%m-%d" if "-" in linha[0] else "%Y/%m/%d"
                    instante_utc = datetime.strptime(linha[0].strip(), formato_data)
                    hora = int("".join(caractere for caractere in linha[1] if caractere.isdigit())[:2])
                    instante_local = (
                        instante_utc.replace(hour=hora, tzinfo=timezone.utc)
                        .astimezone(FUSO)
                        .replace(tzinfo=None)
                    )
                except (ValueError, IndexError):
                    continue
                observacoes[instante_local] = (temperatura, umidade)
    return observacoes


def maximas_diarias(observacoes: dict[datetime, tuple[float, float]]) -> dict[date, float]:
    valores = defaultdict(list)
    for instante, (temperatura, _) in observacoes.items():
        valores[instante.date()].append(temperatura)
    return {dia: max(temperaturas) for dia, temperaturas in valores.items()}


def escrever_csv(caminho: Path, registros: list[dict]) -> None:
    with caminho.open("w", newline="", encoding="utf-8") as destino:
        escrita = csv.DictWriter(destino, fieldnames=list(registros[0]))
        escrita.writeheader()
        for registro in registros:
            escrita.writerow({
                chave: valor.isoformat(sep=" ") if isinstance(valor, datetime)
                else valor.isoformat() if isinstance(valor, date)
                else valor
                for chave, valor in registro.items()
            })


def registros_mmax(
    observacoes: dict[datetime, tuple[float, float]], eventos: list[list[date]], lookup: dict
) -> tuple[list[dict], int]:
    por_dia = {
        dia: (numero, evento[0], evento[-1], len(evento))
        for numero, evento in enumerate(eventos, start=1)
        for dia in evento
    }
    registros, sem_cobertura = [], 0
    for instante in sorted(observacoes):
        if instante.date() not in por_dia or not HORA_INICIO <= instante.hour < HORA_FIM:
            continue
        temperatura, umidade = observacoes[instante]
        chave = (round(temperatura, 1), IAG.arredondar_meio(umidade))
        mmax = lookup.get(chave)
        if mmax is None:
            sem_cobertura += 1
            continue
        numero_evento, inicio, fim, duracao = por_dia[instante.date()]
        registros.append({
            "evento": numero_evento,
            "inicio": inicio,
            "fim": fim,
            "duracao_dias": duracao,
            "data_hora": instante,
            "temperatura_ar_c": temperatura,
            "umidade_relativa_percent": umidade,
            "mmax_met": mmax,
        })
    return registros, sem_cobertura


def frequencias(registros: list[dict]) -> list[dict]:
    contagens = {categoria[2]: 0 for categoria in IAG.CATEGORIAS_ATIVIDADE}
    for registro in registros:
        contagens[IAG.categoria_atividade(registro["mmax_met"])[2]] += 1
    saida = []
    for indice, (met, cor, atividade, rotulo_met) in enumerate(IAG.CATEGORIAS_ATIVIDADE):
        proximo = IAG.CATEGORIAS_ATIVIDADE[indice + 1][0] if indice + 1 < len(IAG.CATEGORIAS_ATIVIDADE) else None
        faixa = f"{rotulo_met} a < {proximo:.1f} MET" if proximo else f"≥ {rotulo_met}"
        saida.append({"atividade": atividade, "faixa": faixa.replace(".", ","), "dias": contagens[atividade], "cor": cor})
    return saida


def tendencia_anual(registros: list[dict]) -> list[dict]:
    valores = defaultdict(list)
    for registro in registros:
        valores[registro["data_hora"].year].append(registro["mmax_met"])
    return [
        {
            "ano": ano,
            "dias": len(serie),
            "mmax_medio_met": sum(serie) / len(serie),
            "mmax_minimo_met": min(serie),
            "mmax_maximo_met": max(serie),
        }
        for ano, serie in sorted(valores.items())
    ]


def main() -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    DADOS_PROCESSADOS.mkdir(exist_ok=True)
    observacoes = ler_observacoes()
    escrever_csv(
        DADOS_PROCESSADOS / "santana-temperatura-umidade-horaria-2006-2025.csv",
        [{"data_hora": instante, "temperatura_ar_c": temperatura, "umidade_relativa_percent": umidade}
         for instante, (temperatura, umidade) in sorted(observacoes.items())],
    )

    maximas = maximas_diarias(observacoes)
    p95 = IAG.percentil(list(maximas.values()), 0.95)
    eventos = IAG.eventos_quentes(maximas, p95)
    eventos_registros = [
        {"evento": indice, "inicio": evento[0], "fim": evento[-1], "duracao_dias": len(evento),
         "tmax_maxima_evento_c": max(maximas[dia] for dia in evento)}
        for indice, evento in enumerate(eventos, start=1)
    ]
    escrever_csv(SAIDA / "santana-eventos-tmax-acima-p95.csv", eventos_registros)

    print(f"Observações válidas: {len(observacoes):,}")
    print(f"P95 da Tmax diária: {p95:.1f} °C")
    print(f"Eventos de pelo menos 3 dias: {len(eventos)}")

    for identificador, (cenario, rotulo_perfil) in CENARIOS.items():
        pasta_perfil = SAIDA / identificador
        pasta_perfil.mkdir(exist_ok=True)
        lookup = IAG.ler_lookup_mmax(ARQUIVO_MODELO, cenario)
        registros, sem_cobertura = registros_mmax(observacoes, eventos, lookup)
        escrever_csv(pasta_perfil / "santana-mmax-horario-por-evento.csv", registros)

        # Versão diária de referência: o menor Mmax das horas diurnas de cada
        # dia quente, equivalente ao gráfico inicial produzido para o IAG.
        por_dia = defaultdict(list)
        for registro in registros:
            por_dia[registro["data_hora"].date()].append(registro)
        diarios = []
        for dia, grupo in sorted(por_dia.items()):
            registro_referencia = grupo[0]
            diarios.append({
                "evento": registro_referencia["evento"],
                "inicio": registro_referencia["inicio"],
                "fim": registro_referencia["fim"],
                "data": dia,
                "duracao_dias": registro_referencia["duracao_dias"],
                "tmax_diaria_c": maximas[dia],
                "mmax_min_met": min(item["mmax_met"] for item in grupo),
                "observacoes_diurnas": len(grupo),
            })
        escrever_csv(pasta_perfil / "santana-mmax-minimo-diario-por-evento.csv", diarios)
        (pasta_perfil / "santana-mmax-minimo-diario-por-evento.svg").write_text(
            IAG.gerar_svg(diarios, p95, min(maximas), max(maximas), estacao="Mirante de Santana",
                           periodo=PERIODO, perfil=rotulo_perfil),
            encoding="utf-8",
        )

        categorias = frequencias(registros)
        with (pasta_perfil / "santana-frequencia-atividades-horaria.csv").open("w", newline="", encoding="utf-8") as destino:
            escrita = csv.DictWriter(destino, fieldnames=["atividade", "faixa", "horas"])
            escrita.writeheader()
            for categoria in categorias:
                escrita.writerow({"atividade": categoria["atividade"], "faixa": categoria["faixa"], "horas": categoria["dias"]})
        (pasta_perfil / "santana-frequencia-atividades-horaria.svg").write_text(
            IAG.gerar_svg_frequencia(categorias, unidade="horas", perfil=rotulo_perfil,
                                     numero_eventos=len(eventos), periodo=PERIODO,
                                     estacao="Mirante de Santana (INMET A701)"),
            encoding="utf-8",
        )

        anual = tendencia_anual(registros)
        inclinacao, intercepto, r2 = IAG.tendencia_linear(anual)
        with (pasta_perfil / "santana-tendencia-mmax-horaria.csv").open("w", newline="", encoding="utf-8") as destino:
            escrita = csv.DictWriter(destino, fieldnames=["ano", "horas", "mmax_medio_met", "mmax_minimo_met", "mmax_maximo_met"])
            escrita.writeheader()
            for linha in anual:
                escrita.writerow({"ano": linha["ano"], "horas": linha["dias"], "mmax_medio_met": linha["mmax_medio_met"],
                                  "mmax_minimo_met": linha["mmax_minimo_met"], "mmax_maximo_met": linha["mmax_maximo_met"]})
        (pasta_perfil / "santana-tendencia-mmax-horaria.svg").write_text(
            IAG.gerar_svg_tendencia(anual, inclinacao, intercepto, r2, unidade="horas", perfil=rotulo_perfil,
                                    ano_inicio=2006, ano_fim=2025, estacao="Mirante de Santana (INMET A701)"),
            encoding="utf-8",
        )
        print(f"{rotulo_perfil}: {len(registros):,} horas avaliadas; {sem_cobertura} sem cobertura da grade; "
              f"tendência {inclinacao * 10:+.3f} MET/década; R²={r2:.3f}")


if __name__ == "__main__":
    main()
