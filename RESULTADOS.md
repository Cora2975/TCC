# Organização do projeto

## Estação IAG

- `dados_estacao_iag/series_temporais_csv/`: séries horárias de entrada.
- `dados_estacao_iag/outputs/`: grade do modelo HEATLim e resultados prévios de comparação do modelo.
- `dados_estacao_iag/resultados/mmax_eventos/adulto_jovem/`: eventos quentes, `Mmax`, frequências de atividade e tendências para adulto jovem.
- `dados_estacao_iag/resultados/mmax_eventos/idosos/`: resultados equivalentes para idosos.
- `dados_estacao_iag/grafico_mmax_eventos_tmax_p95.py`: reproduz os resultados de `Mmax` por evento.

## Mirante de Santana — INMET A701

- `dados_estacao_santana/2006/` a `2025/`: arquivos brutos anuais do INMET.
- `dados_estacao_santana/dados_processados/`: série horária local de temperatura e umidade preparada a partir dos arquivos brutos (`santana-temperatura-umidade-horaria-2006-2025.csv`).
- `dados_estacao_santana/resultados/mmax_eventos/adulto_jovem/`: eventos quentes, `Mmax`, frequências de atividade e tendências para adulto jovem.
- `dados_estacao_santana/resultados/mmax_eventos/idosos/`: resultados equivalentes para idosos.
- `dados_estacao_santana/analise_mmax_santana.py`: reproduz a análise da estação.

## Comparação sazonal entre estações

- `analises_comparativas/mmax_sazonal/`: médias anuais de `Mmax` por estação do ano e respectivas tendências para ambas as estações e faixas etárias.
- `analise_mmax_sazonal.py`: reproduz essa comparação.

Os gráficos estão em SVG; os dados correspondentes estão em CSV na mesma pasta de resultados.
