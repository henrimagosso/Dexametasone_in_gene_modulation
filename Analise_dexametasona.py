import pandas as pd
import os
import numpy as np
from scipy import stats

from Bio import Geo
from Bio import Entrez
from statsmodels.sandbox.stats.multicomp import multipletests


# O arquivo 'GSE52169_series_matrix.txt' deve
# estar na mesma pasta que este script Python
nome_arquivo = 'GSE52169_series_matrix.txt'
separador = '\t'
linhas_para_pular = 62 # Baseado na inspeção manual do cabeçalho

# --- Configuração do Estudo GEO ---
geo_accession = 'GSE52169' # O código GEO do estudo que estamos analisando
# --- CONFIGURAR SEU EMAIL ---
# O NCBI Entrez (usado pelo Biopython para baixar dados) requer um email de contato.
Entrez.email = "henriquemagosso@gmail.com" # <<< SUBSTITUA PELO SEU EMAIL REAL! <<<

# --- Primeira Parte: Lendo o arquivo de dados principais com pandas ---
# Objetivo: Carregar a tabela de expressão (números) em um DataFrame.
print(f"--- Primeira Parte: Carregando dados principais ---")
try:
    print(f"Tentando ler o arquivo: {nome_arquivo}")
    df = pd.read_csv(
        nome_arquivo,
        sep=separador,
        skiprows=linhas_para_pular,
        index_col=0
    )
    df = df.iloc[:-1]

    print("Arquivo de dados principais carregado com sucesso, cabeçalho pulado e rodapé removido!")
    print("\nPrimeiras 5 linhas do DataFrame (df.head()):")
    print(df.head())

    print("\nInformações básicas do DataFrame (df.info()):")
    # Verifique aqui se o índice NÃO inclui '!series_matrix_table_end'
    # e se o número de linhas diminuiu em 1 comparado à leitura anterior sem iloc[:-1]
    df.info() # Print df.info() directly

    print("\nFormato do DataFrame (df.shape):")
    print(df.shape)


except FileNotFoundError:
    print(f"Erro na Primeira Parte: O arquivo '{nome_arquivo}' não foi encontrado.")
    exit()
except Exception as e:
    print(f"Ocorreu um erro na Primeira Parte ao ler o arquivo com pandas ou remover o rodapé: {e}")
    exit()


# --- SEGUNDA PARTE: Extraindo Metadados das Amostras usando Biopython ---
# Objetivo: Usar Biopython para obter metadados de forma robusta do GEO e identificar
# amostras Controle/Dexametasona.

control_samples = [] # Lista para códigos GSM de amostras Controle
dex_samples = []     # Lista para códigos GSM de amostras Dexametasona

print(f"\n--- Segunda Parte: Extraindo Metadados com Biopython ---")
print(f"Conectando ao banco de dados GEO via Entrez para obter metadados do estudo {geo_accession}... (Requer internet)")
print("Isso pode levar alguns segundos...")

try:
    # Busca o registro do estudo (GSE) no banco de dados GEO do NCBI e baixa em formato de texto
    handle = Entrez.efetch(db="geo", id=geo_accession, retmode="text")
    # Usa o módulo Geo do Biopython para parsear (interpretar) o registro baixado
    geo_record = Geo.read(handle) # Geo.read com G minúsculo
    handle.close() # Fecha o handle apos ler
    print(f"Registro GEO {geo_accession} baixado e parseado com sucesso via Biopython.")
    print(f"Número total de amostras encontradas no registro GEO: {len(geo_record.gsm_list)}")

    print("\nClassificando Amostras por Grupo usando Metadados do Biopython...")
    print("Buscando por termos de tratamento ('Vehicle', 'Dex') nas características de cada amostra...")


    df_columns = list(df.columns)

    for gsm in geo_record.gsm_list: # Itera sobre cada objeto Sample (GSM) baixado pelo Biopython
        gsm_code = gsm.accession # O código de acesso da amostra (Ex: 'GSM1254320')

        # --- Lógica de Classificação ---
        # Precisamos procurar nas 'characteristics' por termos que identifiquem Controle ou Dexametasona.
        is_control = False
        is_dex = False

        for char_string in gsm.characteristics:
            cleaned_char = char_string.strip().strip('"').lower()

            # Verificar se é o Controle (buscando pelo termo exato "treatment: vehicle (dmso)")
            if 'treatment: vehicle (dmso)' in cleaned_char:
                is_control = True
                break

            # Verifica se é Dexametasona (buscando por qualquer termo que comece com "treatment: dex")
            # Isso vai capturar "treatment: dex 10 nm", "treatment: dex 100 nm", etc.
            elif cleaned_char.startswith('treatment: dex'):
                 is_dex = True
                 break
            # Nota: Amostras com outros tratamentos (DHT, Enzalutamide, etc.)
            # que não sejam 'Vehicle' nem 'Dex' não serão classificadas aqui

        # Adiciona o código GSM à lista apropriada APENAS SE ele foi classificado
        # E APENAS SE ele for encontrado como uma coluna no DataFrame 'df' (salvaguarda)
        if is_control and gsm_code in df_columns:
            control_samples.append(gsm_code)
        elif is_dex and gsm_code in df_columns:
            dex_samples.append(gsm_code)
        # Amostras com outros tratamentos ou sem tratamento relevante simplesmente não são adicionadas a nenhuma lista


    print("Extração de metadados e classificação das amostras concluída usando Biopython.")

except Exception as e:
    print(f"Ocorreu um erro durante a extração de metadados com Biopython: {e}")
    print("Por favor, verifique se sua conexão com a internet está ativa, se você substituiu 'seu.email@example.com' pelo seu email real, e se instalou Biopython (`pip install biopython`).")
    # Se a extração de metadados falhar, não podemos continuar com a análise diferencial.
    exit()


# --- TERCEIRA PARTE: Verificação Final das Amostras Classificadas ---
# Esta parte usa as listas control_samples e dex_samples populadas pela Segunda Parte (Biopython)
print(f"\n--- Terceira Parte: Resultado da Classificação Final das Amostras ---")
print(f"Códigos GSM para Amostras Controle ({len(control_samples)}):")
print(control_samples)
print(f"\nCódigos GSM para Amostras Dexametasona ({len(dex_samples)}):")
print(dex_samples)


# Verificação básica: Temos amostras nos dois grupos e pelo menos 2 replicatas por grupo (idealmente 3+)?
if not control_samples or not dex_samples:
    print("\nERRO: Não foi possível identificar amostras Controle E Dexametasona suficientes nos metadados via Biopython.")
    print("A análise de Expressão Diferencial não poderá prosseguir sem isso.")
    exit() # Sair se nao encontrar os grupos essenciais
elif len(control_samples) < 2 or len(dex_samples) < 2:
     print("\nAVISO: Idealmente, análise de Expressão Diferencial requer pelo menos 2-3 replicatas por grupo para robustez estatística.")
     print(f"Você identificou {len(control_samples)} amostras Controle e {len(dex_samples)} amostras Dexametasona.")
     print("Vamos continuar com os próximos passos, mas tenha em mente que os resultados podem ser menos confiáveis com poucas replicatas.")
else:
    print("\nAmostras dos grupos Controle e Dexametasona identificadas com sucesso via Biopython.")
    print(f"Identificadas {len(control_samples)} amostras Controle e {len(dex_samples)} amostras Dexametasona.")
    print("Pronto para realizar a Análise de Expressão Diferencial estatística nos dados carregados no DataFrame 'df'.")


# --- QUARTA PARTE: Preparando Dados para Análise Estatística ---
# Objetivo: Selecionar as colunas do DataFrame 'df' que correspondem às amostras Controle e Dexametasona.
print(f"\n--- Quarta Parte: Preparando Dados para Análise Estatística ---")

# Seleciona apenas as colunas do DataFrame que estão nas nossas listas de amostras classificadas
df_control = df[control_samples]
df_dex = df[dex_samples]

print(f"DataFrame com amostras Controle selecionado: {df_control.shape}")
print(f"DataFrame com amostras Dexametasona selecionado: {df_dex.shape}")

# Remover linhas (genes) que tenham valores ausentes (NaN) em qualquer uma das amostras selecionadas
# print("Removendo linhas com valores ausentes...")
# df_control = df_control.dropna()
# df_dex = df_dex.dropna()
# print(f"DataFrames após remover NaNs: Controle {df_control.shape}, Dexametasona {df_dex.shape}")


# --- QUINTA PARTE: Realizar Análise de Expressão Diferencial Estatística ---
# Objetivo: Comparar a expressão de cada gene entre os grupos Controle e Dexametasona.
print(f"\n--- Quinta Parte: Realizando Análise de Expressão Diferencial (Teste T) ---")

p_values = []
t_statistics = []
mean_control_values = []
mean_dex_values = []
genes = [] # Para guardar o ID do gene/sonda na ordem dos resultados

# axis=1 para calcular a média/aplicar teste ao longo das colunas (para cada linha/gene)
for index, row in df.iterrows():
    gene = index # O ID do gene/sonda é o índice da linha
    genes.append(gene)

    # Obtém os valores de expressão para esta linha (gene) nos grupos Controle e Dexametasona
    # Use .values para obter um array numpy, e dropna() para remover NaNs para o teste
    control_values = row[control_samples].dropna().values
    dex_values = row[dex_samples].dropna().values

    # Calcula as médias para o Fold Change depois
    mean_control = np.mean(control_values) if len(control_values) > 0 else 0
    mean_dex = np.mean(dex_values) if len(dex_values) > 0 else 0
    mean_control_values.append(mean_control)
    mean_dex_values.append(mean_dex)

    # 'nan_policy="omit"' ignora NaNs nos arrays passados para o teste
    # 'equal_var=False' realiza o teste t de Welch, que não assume variâncias iguais (geralmente mais apropriado para dados biológicos)
    if len(control_values) >= 2 and len(dex_values) >= 2 and np.std(control_values) > 1e-9 and np.std(dex_values) > 1e-9: # Garante min 2 dados e variancia > 0
        try:
            t_stat, p_val = stats.ttest_ind(
                control_values,
                dex_values,
                equal_var=False, # Welch's t-test
                nan_policy='omit' # Omite NaNs específicos da scipy
            )
            t_statistics.append(t_stat)
            p_values.append(p_val)
        except Exception as e:
             # Em casos raros, o teste t pode falhar mesmo com dados suficientes (ex: todos os valores iguais apos dropna)
             print(f"Aviso: Teste t falhou para o gene {gene}: {e}. Atribuindo NaN para estatísticas.")
             t_statistics.append(np.nan)
             p_values.append(np.nan)
    else:
        t_statistics.append(np.nan)
        p_values.append(np.nan)

print("Análise de Expressão Diferencial (Testes T) concluída.")

# --- SEXTA PARTE: Corrigir Múltiplos Testes e Calcular Fold Change ---
# Objetivo: Ajustar os p-valores e calcular a magnitude da mudança na expressão.
print(f"\n--- Sexta Parte: Corrigindo Múltiplos Testes e Calculando Fold Change ---")

# Corrigir p-valores para múltiplos testes (FDR)
# Ignoramos p-valores que são NaN no calculo da correcao
p_values_valid = np.array(p_values)
mask_valid_p = ~np.isnan(p_values_valid) # Mascara para valores nao ausentes
p_values_valid_filtered = p_values_valid[mask_valid_p]

if len(p_values_valid_filtered) > 0:
    # results é uma tupla: (reject, pvals_corrected, _, _)
    reject, pvals_corrected, _, _ = multipletests(p_values_valid_filtered, method='fdr_bh')

    # Criar um array de p-valores ajustados com NaNs onde os p-valores originais eram NaN
    adjusted_p_values = np.full(len(p_values), np.nan)
    adjusted_p_values[mask_valid_p] = pvals_corrected
else:
     # Se nao houver p-valores validos, todos os ajustados sao NaN
     adjusted_p_values = np.full(len(p_values), np.nan)

fold_changes = (np.array(mean_dex_values) + 1) / (np.array(mean_control_values) + 1)
log2_fold_changes = np.log2(fold_changes) # log base 2

print("Correção de múltiplos testes e cálculo de Fold Change concluídos.")


# --- SÉTIMA PARTE: Compilar Resultados e Identificar Genes Significativos ---
# Objetivo: Reunir todas as estatísticas em um DataFrame de resultados e filtrar.
print(f"\n--- Sétima Parte: Compilando Resultados e Identificando Genes Significativos ---")

# Cria um DataFrame com todos os resultados
results_df = pd.DataFrame({
    'Mean_Control': mean_control_values,
    'Mean_Dexamethasone': mean_dex_values,
    'Fold_Change': fold_changes,
    'Log2_Fold_Change': log2_fold_changes,
    'T_Statistic': t_statistics,
    'P_Value': p_values,
    'Adjusted_P_Value (FDR)': adjusted_p_values
}, index=genes) # Usa os IDs dos genes/sondas como índice

# Ordena os resultados pelo p-valor ajustado
results_df = results_df.sort_values(by='Adjusted_P_Value (FDR)')

print("\nResultados completos (primeiras 10 linhas ordenadas por p-valor ajustado):")
print(results_df.head(10))

# --- Filtrar por Significância ---
# Define os limites para significância estatística e mudança biológica
# Limite comum: p-valor ajustado < 0.05 E |Log2(Fold Change)| >= 1 (ou seja, Fold Change >= 2 ou <= 0.5)
fdr_threshold = 0.05
log2fc_threshold = 1.0 # Corresponde a um Fold Change de 2x

# Filtra o DataFrame de resultados
significant_genes_df = results_df[
    (results_df['Adjusted_P_Value (FDR)'] < fdr_threshold) &
    (np.abs(results_df['Log2_Fold_Change']) >= log2fc_threshold)
]

print(f"\nIdentificados {len(significant_genes_df)} genes/sondas significativamente expressos:")
print(f"(Limite: p-valor ajustado < {fdr_threshold} e |Log2(Fold Change)| >= {log2fc_threshold})")

# Exibe as primeiras 10 linhas dos genes significativos encontrados
print("\nGenes/Sondas Significativos (primeiras 10 linhas):")
print(significant_genes_df.head(10))


# --- OITAVA PARTE: Visualização de Resultados (Próximos Passos) ---
# O código para gerar Volcano Plot, Heatmap, etc., virá depois.
# Geralmente, para visualização, precisaríamos instalar matplotlib e/ou seaborn.
# import matplotlib.pyplot as plt # Se for usar matplotlib
# import seaborn as sns # Se for usar seaborn

print("\n--- Oitava Parte: Visualização de Resultados (Próximos Passos) ---")
print("Para visualizar (Volcano Plot, Heatmap), instale matplotlib e seaborn (`pip install matplotlib seaborn`).")
print("O código para visualização virá aqui.")

# Exemplo BÁSICO de como seria uma visualização (apenas esboço, requer matplotlib)
# if not significant_genes_df.empty:
#    plt.figure(figsize=(10, 6))
#    sns.scatterplot(data=results_df, x='Log2_Fold_Change', y=-np.log10(results_df['Adjusted_P_Value (FDR)']), alpha=0.5, label='Não Significativo')
#    sns.scatterplot(data=significant_genes_df, x='Log2_Fold_Change', y=-np.log10(significant_genes_df['Adjusted_P_Value (FDR)']), color='red', alpha=0.6, label='Significativo')
#    plt.axhline(-np.log10(fdr_threshold), color='grey', linestyle='--', label=f'Limite FDR ({fdr_threshold})')
#    plt.axvline(log2fc_threshold, color='grey', linestyle='--', label=f'Limite Log2FC ({log2fc_threshold})')
#    plt.axvline(-log2fc_threshold, color='grey', linestyle='--')
#    plt.title('Volcano Plot da Análise de Expressão Diferencial')
#    plt.xlabel('Log2(Fold Change)')
#    plt.ylabel('-Log10(p-valor Ajustado)')
#    plt.legend()
#    plt.grid(True, linestyle='--', alpha=0.6)
#    plt.show()


print("\nScript concluído.")
