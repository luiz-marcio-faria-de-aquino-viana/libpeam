"""
COPPE/UFRJ 29/SET/2020)
Probabilidade e Estatistica para Aprendizado de Maquina
 
Lista 4 - Lista de Exercicios
 
Nome: Luiz Marcio Faria de Aquino Viana
DRE: 120048833
CPF: 024.723.347-10
RG: 08855128-8 IFP-RJ

Professores:
    Rosa Leao
    Edmundo de Souza e Silva
    Daniel Menasche
 
LIBPEAM: Biblioteca de Probabilidade e Estatistica para Aprendizado de Maquina
"""
import sys as s
import math as M
import statistics as S
import matplotlib.pyplot as plt

#from scipy.stats import geom as F_geom
#from scipy.stats import binom as F_binom
#from scipy.stats import poisson as F_poisson
#from scipy.stats import expon as F_exp
#from scipy.stats import norm as F_normal
#from scipy.stats import weibull_min as F_weibull
from scipy.stats import lognorm as F_lognorm
from scipy.stats import beta as F_beta
from scipy.stats import chi2 as F_chi2
from scipy.stats import f as F_f

# peam_readdata(): obtem uma lista a partir de um arquivo de dados (sem cabecalho)
# fileName - nome do arquivo de dados (CSV: Xn,Yn)
# retorno: lista com os elementos do arquivo
def peam_readdata(fileName):
    # VARS DEF
    #
    f = None                    # objeto de acesso ao arquivo
    N = 0                       # numero de linhas no arquivo
    buf_data = ''               # buffer de entrada de dados
    buf_data = ''               # buffer de dados de entrada (arquivo completo)
    tmp_ls = []                 # lista de linhas de dados de entrada (\n)
    tmp_buf = ''                # 1 linha de dados de entrada (,)
    tmp_data = []               # lista de dados de 1 linha de entrada
    X = ''                      # valor X da entrada (string) 
    Y = 0.0                     # valor Y da entrada (double) 
    dst_data = []               # lista com os dados de 1 elemento [Xn,Yn]
    dst_ls = []                 # retorno: lista de elementos do arquivo (sem cabecalho)
    
    # EXEC
    print("ARQUIVO DE ENTRADA: ", fileName, '\n')

    f = open(fileName, 'r')
    buf_data = f.read(65536)
    f.close()

    tmp_ls = buf_data.split('\n')
    N = 0
    for tmp_buf in tmp_ls:
        tmp_data = tmp_buf.split(',')
        
        if(N == 0):
            print("COL[0]: ", tmp_data[0])
            print("COL[1]: ", tmp_data[1])
        else:
            #print("COL[0]: ", tmp_data[0])
            #print("COL[1]: ", tmp_data[1])
            X = tmp_data[0]
            Y = float(tmp_data[1])
            dst_data = []
            dst_data.append(X)
            dst_data.append(Y)
            dst_ls.append(dst_data)
        N = N + 1            

    # RET
    return dst_ls


# peam_qtd(): calcula a quantidade de elementos em uma lista
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# retorno: numero de elementos
def peam_qtd(src_ls):
    # VARS DEF
    #
    src_data = ['',0.0]         # lista com os dados de 1 elemento [Xn,Yn]
    N = 0.0                     # retorno: numero de elementos
    
    # EXEC
    for src_data in src_ls:
        N = N + 1.0

    # RET
    return N    

# peam_soma(): obtem a soma dos valores da lista
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# retorno: soma dos elementos da lista
def peam_soma(src_ls):
    # VARS DEF
    #
    src_data = ['',0.0]         # lista com os dados de 1 elemento [Xn,Yn]
    soma = 0.0                  # retorno: soma dos elementos da lista
    
    # EXEC
    for src_data in src_ls:
        soma = soma + src_data[1]

    # RET
    return soma

# peam_min(): obtem o menor valor da lista
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# retorno: menor elemento da lista
def peam_min(src_ls):
    # VARS DEF
    #
    src_data = ['',0.0]         # lista com os dados de 1 elemento [Xn,Yn]
    min = 0.0                   # retorno: menor elemento da lista
    
    # EXEC
    src_data = src_ls[0]
    min = src_data[1]
    
    for src_data in src_ls:
        if(src_data[1] < min):
            min = src_data[1]

    # RET
    return min

# peam_max(): obtem o maior valor da lista
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# retorno: maior elemento da lista
def peam_max(src_ls):
    # VARS DEF
    #
    src_data = ['',0.0]         # lista com os dados de 1 elemento [Xn,Yn]
    min = 0.0                   # retorno: maior elemento da lista
    
    # EXEC
    src_data = src_ls[0]
    max = src_data[1]
    
    for src_data in src_ls:
        if(src_data[1] > max):
            max = src_data[1]

    # RET
    return max

# peam_media(): calcula a media dos elementos em uma lista
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# retorno: valor calculado da media
def peam_media(src_ls):
    # VARS DEF
    #
    src_data = ['',0.0]         # lista com os dados de 1 elemento [Xn,Yn]
    N = 0.0                     # numero de elementos
    soma = 0.0                  # somatorio dos elementos
    media = 0.0                 # retorno: valor calculado da media
    
    # EXEC
    for src_data in src_ls:
        #print("X: ", src_data[0])
        #print("Y: ", src_data[1])
        soma = soma + src_data[1]
        N = N + 1.0
    media = soma / N       

    # RET
    return media    
    
# peam_var(): calcula a variancia dos elementos em uma lista
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# media - media dos elementos da lista
# retorno: valor calculado da variancia
def peam_var(src_ls, media):
    # VARS DEF
    #
    src_data = [0.0,0.0]        # lista com os dados de 1 elemento [Xn,Yn]
    N = 0.0                     # numero de elementos
    dist_media2 = 0.0           # quandrado da distancia a media
    soma2 = 0.0                 # somatorio do quandrado da distancia a media
    var = 0.0                   # retorno: valor calculado da variancia
    
    # EXEC
    for src_data in src_ls:
        dist_media2 = (src_data[1] - media) * (src_data[1] - media)
        soma2 = soma2 + dist_media2 
        N = N + 1.0
    var = soma2 / N       

    # RET
    return var    
    
# peam_std(): calcula o desvio padrao dos elementos em uma lista
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# var - variancia dos elementos da lista
# retorno: valor calculado do desvio padrao
def peam_std(var):
    # VARS DEF
    #
    std = 0.0                   # retorno: valor calculado do desvio padrao
    
    # EXEC
    std = M.sqrt(var)       

    # RET
    return std    
    
# peam_outliers(): lista os elementos que estao alem de Z desvios padroes da media
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# std - desvio padrao dos elementos da lista
# Z - numero de desvios padroes em relacao a media
# bNegativo - considera valores negativos (True=SIM / False=NAO)
# retorno: lista dos elementos distantes Z desvios padores da media 
def peam_outliers(src_ls, media, std, Z, bNegativo):
    # VARS DEF
    #
    src_data = [0.0,0.0]        # lista com os dados de 1 elemento [Xn,Yn]
    Z_data = 0.0                # distancia do elemento a media
    val = 0.0                   # valor do dado
    N_in = 0.0                  # numero de itens no workset
    N_out = 0.0                 # numero de outliers
    dst_ls = []                 # retorno: lista dos outliers
    
    # EXEC
    for src_data in src_ls:
        val = src_data[1]
        if(bNegativo == True):
            #print("ACEITA_NEGATIVO")
            Z_data = abs((val - media) / std)
            if(Z_data > Z):
                dst_ls.append(src_data)
                N_out = N_out + 1
            else:
                N_in = N_in + 1
        else:
            #print("NAO_ACEITA_NEGATIVO")
            if(val >= 0):
                #print("POSITIVO", val)
                Z_data = abs((val - media) / std)
                if(Z_data > Z):
                    dst_ls.append(src_data)
                    N_out = N_out + 1
                else:
                    N_in = N_in + 1
            else:
                #print("NEGATIVO", val)
                dst_ls.append(src_data)
                N_out = N_out + 1

    #print("WORKSET_SIZE:", N_in)
    #print("OUTLIERS_SIZE:", N_out)
        
    # RET
    return dst_ls    
    
# peam_workset(): lista os elementos que estao abaixo de Z desvios padroes da media
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# std - desvio padrao dos elementos da lista
# Z - numero de desvios padroes em relacao a media
# bNegativo - considera valores negativos (True=SIM / False=NAO)
# retorno: lista dos elementos distantes Z desvios padores da media 
def peam_workset(src_ls, media, std, Z, bNegativo):
    # VARS DEF
    #
    N_in = 0.0                  # numero de itens no workset
    N_out = 0.0                 # numero de outliers
    src_data = [0.0,0.0]        # lista com os dados de 1 elemento [Xn,Yn]
    Z_data = 0.0                # distancia do elemento a media
    dst_ls = []                 # retorno: lista dos outliers
    
    # EXEC
    for src_data in src_ls:
        if(bNegativo == True):
            Z_data = abs((src_data[1] - media) / std)
            if(Z_data <= Z):
                dst_ls.append(src_data)
                N_in = N_in + 1
            else:
                N_out = N_out + 1
        else:
            if(src_data[1] >= 0):
                Z_data = abs((src_data[1] - media) / std)
                if(Z_data <= Z):
                    dst_ls.append(src_data)
                    N_in = N_in + 1
                else:
                    N_out = N_out + 1
            else:
                N_out = N_out + 1
                
    #print("WORKSET_SIZE:", N_in)
    #print("OUTLIERS_SIZE:", N_out)
        
    # RET
    return dst_ls    
    
# peam_adiciona_freq(): obtem a faixa de ocorrencia pelo valor
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# val - valor a ser adicionado a frequencia de valores da faixa
# N - numero de elementos de entrada
# retorno: lista atualizada de faixas 
def peam_adiciona_freq(faixa_ls, val, N):
    for faixa_data in faixa_ls:
        if( (val >= faixa_data[1]) and (val <= faixa_data[2]) ):
            faixa_data[3] = faixa_data[3] + 1
            faixa_data[4] = faixa_data[3] / N
    return faixa_ls        
    
# peam_lista_freq(): lista as frequencias de ocorrencia por faixa de valores
# src_ls - lista de elementos [ [X1,Y1],[X2,Y2], ... ,[Xn,Yn] ]
# n - numero de faixas de valores
# min - valor minimo da lista
# max - valor maximo da lista
# N - numero de elementos de entrada
def peam_lista_freq(src_ls, n, min, max, N):
    # VARS DEF
    #
    src_data = []                   # elemento da lista de entrada
    faixa_data = []                 # item da lista de faixas
    faixa_val = (max - min) / n     # tamanho da faixa
    faixa_freq = 0.0                # frequencia inicial das faixas
    faixa_prob = 0.0                # probabilidade de ocorrencia de cada faixa
    lim_min = 0.0                   # limite inferior da faixa
    lim_max = 0.0                   # limite superior da faixa
    val_curr = 0.0                  # valor da faixa atual
    faixa_ls = []                   # retorno: lista de faixas
    
    # EXEC
    val_curr = 0.0
    for i in range(int(n)):
        lim_min = val_curr
        lim_max = val_curr + faixa_val

        faixa_data = []
        faixa_data.append(i)
        faixa_data.append(lim_min)
        faixa_data.append(lim_max)
        faixa_data.append(faixa_freq)
        faixa_data.append(faixa_prob)

        faixa_ls.append(faixa_data)

        val_curr = val_curr + faixa_val

    for src_data in src_ls:
        faixa_ls = peam_adiciona_freq(faixa_ls, src_data[1], N) 
                
    # RET
    return faixa_ls    

# find_item(): funcao que pesquisa um item da lista
# key - chave de classificacao do item
# ls - lista de elementos para pesquisa
# retorno: elemento da lista contendo a chave correspondente
def find_item(key, ls):
    # VARS DEF
    #
    data_ls = []
    
    # EXEC
    #
    for data_ls in ls:
        if data_ls[0] == key:
            return data_ls

    # RET
    return []

# uni_workset(workset1, workset2)
# workset_ls1 - lista de dados da primeira amostra
# workset_ls2 - lista de dados da segunda amostra
# retorno: lista com os dois worksets combinados [key,val]
def uni_workset(workset_ls1, workset_ls2):
    # VARS DEF
    #
    key1 = ''
    key2 = ''
    val1 = 0.0
    val2 = 0.0
    val = 0.0
    workset_data1 = []
    workset_data2 = []
    result_data = []
    result_ls = []

    # EXEC
    #
    for workset_data1 in workset_ls1:
        key1 = workset_data1[0]
        val1 = workset_data1[1]
        
        workset_data2 = find_item(key1, workset_ls2)
        if len(workset_data2) > 0:
            val2 = workset_data2[1]
        else:
            val2 = 0.0
            
        val = val1 + val2

        result_data = []
        result_data.append(key1)
        result_data.append(val)
        
        result_ls.append(result_data)
        
    for workset_data2 in workset_ls2:
        key2 = workset_data2[0]
        val2 = workset_data2[1]
        
        workset_data1 = find_item(key2, workset_ls1)
        if len(workset_data1) == 0:
            result_data = []
            result_data.append(key2)
            result_data.append(val2)
        
            result_ls.append(result_data)

    # RET
    return result_ls        

# add_xval_workset(workset1, workset2)
# workset_ls1 - lista de dados da primeira amostra
# workset_lsC - lista de dados da segunda amostra completa
# retorno: lista com os dados do primeiro workset com as chaves que faltavam =[key,val]
def add_xval_workset(workset_ls1, workset_lsC):
    # VARS DEF
    #
    keyC = ''
    workset_dataC = []
    workset_data1 = []
    result_data = []
    result_ls = []

    for workset_data1 in workset_ls1:
        result_ls.append(workset_data1)

    for workset_dataC in workset_lsC:
        keyC = workset_dataC[0]
        
        workset_data1 = find_item(keyC, workset_ls1)
        if len(workset_data1) == 0:
            result_data = []
            result_data.append(keyC)
            result_data.append(0.0)
        
            result_ls.append(result_data)

    # RET
    return result_ls        

# add_valgeom_faixa(faixa_ls1)
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# retorno: lista de faixas com a probabilidade obtida da Disrtibuicao Geometrica
#       tfaixa_ls - [ [faixa1,min1,max1,freq1,p1,f1,E1],[faixa2,min2,max2,freq2,p2,f2,E2], ... ,[faixaN,minN,maxN,freqN,pN,fN,EN] ]
def add_valgeom_faixa(faixa_ls1, C1, pS):
    # VARS DEF
    #
    faixa_data1 = []
    faixa_meio1 = 0.0
    fx = 0.0
    result_data = []
    result_ls = []    
    freqN = 0.0
    fx = 0.0
    Ei = 0.0

    for faixa_data1 in faixa_ls1:
        freqN = freqN + faixa_data1[3]    
    
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        fx = peam_func_fpm_geom(faixa_meio1, C1, pS)
        Ei = fx * freqN 

        result_data = []
        result_data.append(faixa_data1[0])
        result_data.append(faixa_data1[1])
        result_data.append(faixa_data1[2])
        result_data.append(faixa_data1[3])
        result_data.append(faixa_data1[4])
        result_data.append(fx)
        result_data.append(Ei)
        result_data.append(freqN)
        
        result_ls.append(result_data)

    # RET
    return result_ls     

def hist_interval(n, std, min1, max1):
    # VARS DEF
    m = M.floor(1.0 + 3.3 * M.log10(n))   # numero de intervalos
    sz = (min1 - max1) / m                # tamanho do intervalo na amostra
    p_sz = 1.0 / m                        # tamanho do intervalo regularizado [0.1]
    exp = - 1.0 / 3.0                     # expoente para o calculo do "bin"
    b = 3.49 * std * M.pow(n, exp)        # valor calculado do "bin"
    
    result_ls = []
    result_ls.append(m)
    result_ls.append(b)
    result_ls.append(sz)
    result_ls.append(p_sz)

    # RET
    return result_ls

# FUNCOES ESTATISTICAS
#

# peam_func_erf(): calcula a funcao erro de Gauss
# x - valor de entrada
# retorno: valor de retorno da funcao erro de Gauss
def peam_func_erf(t):
    # VARS DEF
    #
    erf = 0.0
    
    # EXEC
    #        
    erf = - 2 * t

    # RET
    return erf
    
# peam_func_fpm_geom(): calcula a probabilidade de ocorrencia de um valor em uma distribuicao geometrica
# k - numero de falhas ate a ocorrencia de sucesso
# C1 - constante de ajuste
# pS - probabilidade de sucesso
# retorno: probabilidade P[X=k|p]
def peam_func_fpm_geom(k, C1, pS):
    # VARS DEF
    #
    pF = (1.0 - pS)
    Px = C1 * M.pow(pF, k - 1) * pS

    # RET
    return Px

# peam_func_fpa_geom(): calcula a probabilidade P[X<x] em uma distribuicao geometrica
# k - numero de falhas ate a ocorrencia de sucesso
# C1 - constante de ajuste
# pS - probabilidade de sucesso
# retorno: probabilidade P[X<=k|p]
def peam_func_fpa_geom(k, C1, pS):
    # VARS DEF
    #
    PxAcum = 0.0
    
    x = 0.0
    while x <= k:
        Px = peam_func_fpm_geom(x, C1, pS)      
        PxAcum = PxAcum + Px
        x = x + 1.0
            
    # RET
    return PxAcum

# FUNCOES DE DISTRIBUICAO
#

# F_weibull
#

# peam_func_pdf_weibull(): calcula a probabilidade de ocorrencia de um valor em uma distribuicao F_weibull
# x - variavel aleatoria
# C1 - constante de ajuste
# val_lambda - parametro "lambda" da distribuicao
# k - parametro "k" (forma da curva; k<1 = exponencial acentuada; k=1 = exponencial; k>1 = polinomial acentuada) 
# retorno: probabilidade P[X=x|k]
def peam_func_pdf_weibull(x, C1, val_lambda, k):
    # EXEC
    #
    p1 = k / val_lambda
    
    p2 = x / val_lambda
    
    p3 = - M.pow(p2, k)
    
    fx = C1 * p1 * M.pow(p2, k - 1) * M.exp(p3)

    # RET
    return fx

# peam_func_cdf_weibull(): calcula a probabilidade acumulada de ocorrencia de um valor em uma distribuicao F_weibull
# x - variavel aleatoria
# C1 - constante de ajuste
# val_lambda - parametro "lambda" da distribuicao
# k - parametro "k" (forma da curva; k<1 = exponencial acentuada; k=1 = exponencial; k>1 = polinomial acentuada) 
# retorno: probabilidade P[X<x|k]
def peam_func_cdf_weibull(x, C1, val_lambda, k):
    # EXEC
    #
    p1 = x / val_lambda
    
    p2 = - M.pow(p1, k)
    
    Fx = C1 * (1 - M.exp(p2))

    # RET
    return Fx

# F_lognorm
#

# peam_func_pdf_lognorm(): calcula a probabilidade de ocorrencia de um valor em uma distribuicao F_lognorm
# x - variavel aleatoria
# C1 - constante de ajuste
# media - media da distribuicao
# std - desvio padrao da distribuicao 
# retorno: probabilidade P[X=x|media,std]
def peam_func_pdf_lognorm(x, C1, media, std):
    # EXEC
    #
    #print(media)
    scale = M.exp(media)        # escala da curva
    loc = 0.0                   # deslocamento em "x"
    fx = C1 * F_lognorm.pdf(x, std, loc, scale)

    # RET
    return fx

# peam_func_cdf_lognorm(): calcula a probabilidade acumulada de ocorrencia de um valor em uma distribuicao F_lognorm
# x - variavel aleatoria
# C1 - constante de ajuste
# media - media da distribuicao
# std - desvio padrao da distribuicao 
# retorno: probabilidade acumulada P[X<=x|a,b]
def peam_func_cdf_lognorm(x, C1, media, std):
    # EXEC
    #
    #print(media)
    scale = M.exp(media)        # escala da curva
    loc = 0.0                   # deslocamento em "x"
    fx = C1 * F_lognorm.cdf(x, std, loc, scale)

    # RET
    return fx

# F_beta
#

# peam_func_pdf_beta(): calcula a probabilidade de ocorrencia de um valor em uma distribuicao F_beta
# x - variavel aleatoria
# C1 - constante de ajuste
# a - parametro "a" da distribuicao 
# b - parametro "b" da distribuicao
# retorno: probabilidade P[X=x|a,b]
def peam_func_pdf_beta(x, C1, a, b):
    # EXEC
    #
    Px = C1 * F_beta.pdf(x, a, b)

    # RET
    return Px

# peam_func_cdf_beta(): calcula a probabilidade acumulada de ocorrencia de um valor em uma distribuicao F_beta
# x - variavel aleatoria
# C1 - constante de ajuste
# a - parametro "a" da distribuicao 
# b - parametro "b" da distribuicao
# retorno: probabilidade acumulada P[X<=x|a,b]
def peam_func_cdf_beta(x, C1, a, b):
    # EXEC
    #
    Px = C1 * F_beta.cdf(x, a, b)

    # RET
    return Px

# F_exp
#

# peam_func_pdf_exp(): calcula a probabilidade de ocorrencia de um valor em uma distribuicao exponencial
# x - variavel aleatoria
# C1 - constante de ajuste
# lambda - valor do parametro lambda da distribuicao
# retorno: valor da funcao de distribuicao de probabilidade para f(X=x|lambda)
def peam_func_pdf_exp(x, C1, val_lambda):
    # EXEC
    #
    fx = C1 * (val_lambda * x) / M.exp(val_lambda * x)

    # RET
    return fx

# peam_func_cdf_exp(): calcula a probabilidade acumulada de um valor em uma distribuicao exponencial
# x - variavel aleatoria
# C1 - constante de ajuste
# lambda - valor do parametro lambda da distribuicao
# retorno: probabilidade acumulada F(X<=x|lambda)
def peam_func_cdf_exp(x, C1, val_lambda):
    # EXEC
    #
    Fx = 1.0 - (1.0 / M.exp(val_lambda * x))

    # RET
    return Fx

# F_norm
#

# peam_func_pdf_norm(): calcula a probabilidade de ocorrencia de um valor em uma distribuicao normal
# x - variavel aleatoria
# C1 - constante de ajuste
# media - media da distribuicao
# std - desvio padrao da distribuicao
# retorno: valor da funcao de distribuicao para f(X=x|media,std)
def peam_func_pdf_norm(x, C1, media, std):
    # EXEC
    #
    f = S.NormalDist(media, std)
    fx = C1 * f.pdf(x)
    
    # RET
    return fx

# peam_func_cdf_norm(): calcula a probabilidade acumulada de ocorrencia de um valor em uma distribuicao normal
# x - variavel aleatoria
# C1 - coeficiente de ajuste
# media - media da distribuicao
# std - desvio padrao da distribuicao
# retorno: probabilidade f[X<x|media,std]
def peam_func_cdf_norm(x, C1, media, std):
    # EXEC
    #
    f = S.NormalDist(media, std)
    Fx = C1 * f.cdf(x)
    
    # RET
    return Fx

# F_f
#

# peam_func_pdf_F(): calcula a probabilidade de ocorrencia de um valor em uma distribuicao F (Fisher-Snedecor)
# x - variavel aleatoria
# d1 - parametro "d1" da distribuicao 
# d2 - parametro "d2" da distribuicao
# retorno: probabilidade P[X=x|a1,b1,N,p]
def peam_func_pdf_F(x, d1, d2):
    # VARS DEF
    #
    f = F_f(d1, d2)
    fx = f.pdf(x)

    # RET
    return fx

# peam_func_cdf_F(): calcula a probabilidade acumulada de um valor em uma distribuicao F
# x - variavel aleatoria
# d1 - valor do parametro "d1" da distribuicao
# d2 - valor do parametro "d2" da distribuicao
# retorno: probabilidade f[X<x|d1,d2]
def peam_func_cdf_F(x, d1, d2):
    # VARS DEF
    #
    f = F_f(d1, d2)
    Fx = f.cdf(x)

    # RET
    return Fx

# F_chi2
#

# peam_func_pdf_chi2(): calcula a probabilidade de ocorrencia de um valor em uma distribuicao Chi-Square
# x - variavel aleatoria
# df - parametro "df" da distribuicao 
# retorno: probabilidade P[X=x|df]
def peam_func_pdf_chi2(x, df):
    # VARS DEF
    #
    f = F_chi2(df)
    fx = f.pdf(x)

    # RET
    return fx

# peam_func_cdf_chi2(): calcula a probabilidade acumulada de um valor em uma distribuicao Chi-Square
# x - variavel aleatoria
# df - valor do parametro "df" da distribuicao
# retorno: probabilidade f[X<x|df]
def peam_func_cdf_chi2(x, df):
    # VARS DEF
    #
    f = F_chi2(df)
    Fx = f.cdf(x)

    # RET
    return Fx

# CALCULO DO VALOR DE "X" EM FUNCAO DA AREA
#

# peam_valx_chi2(): calcula o valor de "X" em funcao da area sobre a curva
# area - valor da area sobre a curva
# df - valor do parametro "df" da distribuicao
# Niter - numero maximo de iteracoes
# retorno: valor aproximado de "X"
def peam_valx_chi2(area, df, min_x, max_x):
    # VARS DEF
    #
    result_ls = []
    
    Niter = 100000.0
    
    dif = max_x - min_x    
    
    dx = dif / Niter
    
    x_pos = 0.5
    
    x = x_pos * dx
    
    Fx = peam_func_cdf_chi2(x, df) 
    
    while (x_pos < Niter) and (Fx < area):
        x = x_pos * dx
        Fx = peam_func_cdf_chi2(x, df)         
        
        #print("x_pos:", x_pos, "x:", x, "Fx:", Fx)

        x_pos = x_pos + 1.0

    area_dif = abs(Fx - area)        
        
    result_ls.append(x)
    result_ls.append(area)
    result_ls.append(Fx)
    result_ls.append(area_dif)

    # RET
    return result_ls

# peam_valx_F(): calcula o valor de "X" em funcao da area sobre a curva
# area - valor da area sobre a curva
# d1 - valor do parametro "d1" da distribuicao
# d2 - valor do parametro "d2" da distribuicao
# Niter - numero maximo de iteracoes
# retorno: valor aproximado de "X"
def peam_valx_F(area, d1, d2, min_x, max_x):
    # VARS DEF
    #
    result_ls = []
    
    Niter = 100000.0
    
    dif = max_x - min_x    
    
    dx = dif / Niter
    
    x_pos = 0.5
    
    x = x_pos * dx
    
    Fx = peam_func_cdf_F(x, d1, d2) 
    
    while (x_pos < Niter) and (Fx < area):
        x = x_pos * dx
        Fx = peam_func_cdf_F(x, d1, d2)         
        
        #print("x_pos:", x_pos, "x:", x, "Fx:", Fx)

        x_pos = x_pos + 1.0

    area_dif = abs(Fx - area)        
        
    result_ls.append(x)
    result_ls.append(area)
    result_ls.append(Fx)
    result_ls.append(area_dif)

    # RET
    return result_ls

# MLE
#

# peam_mle_geom(): calcula o valor do parametro "p" em funcao dos valores das frequencias das faixas
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# retorno: lista contendo a constante de ajuste e a probabilidade calculadas [C1,pS,pF]
def peam_mle_geom(faixa_ls):
    # VARS DEF
    #
    faixa_data = []                 # elemento da lista de faixas [Nfaixa,min,max,freq]
    faixa_meio = 0.0                # valor central da faixa
    faixa_prob = 0.0                # valor central da faixa
    valor_exp = 0.0                 # valor experado    
    pS = 0.0                        # parametro calculado - pS (sucesso)
    pF = 0.0                        # parametro calculado - pF (falha)
    C1 = 0.0                        # constante calculada - C1
    fx = 0.0                        # valor da funcao distribuicao
    somaF = 0.0                     # somatorio dos valores da funcao de distribuicao
    result_ls = []                  # resultado: lista com valores calculados para "C1" e "p"
    
    # EXEC
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_prob = faixa_data[4]; 
        valor_exp = valor_exp + (faixa_meio * faixa_prob)

    # Distribuicao Geometrica 
    # E(X) = 1/P(X)
    # VAR(X) = (1 - P(X)) / (P(X) * P(X))
    
    pS = 1.0 / valor_exp;
    pF = 1.0 - pS

    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        fx = peam_func_fpm_geom(faixa_meio, 1.0, pS)
        somaF = somaF + fx
    
    C1 = 1.0 / somaF
    
    result_ls.append(C1)
    result_ls.append(pS)
    result_ls.append(pF)
    
    # RET
    return result_ls

# peam_mle_beta(): calcula o valor do parametro "a" e "b" em funcao dos valores das frequencias das faixas
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# retorno: lista com valores calculados dos parametros [C1,a, b]
def peam_mle_beta(faixa_ls):
    # VARS DEF
    #
    faixa_data = []                 # elemento da lista de faixas [Nfaixa,min,max,freq,p]
    faixa_meio = 0.0                # valor central da faixa
    faixa_prob = 0.0                # valor central da faixa
    valor_exp = 0.0                 # valor experado    
    valor_exp2 = 0.0                # valor experado de E[X2]    
    valor_var = 0.0                 # valor da variancia
    valor_v = 0.0                   # valor calculado - v
    C1 = 0.0                        # parametro calculado - C1
    a = 0.0                         # parametro calculado - alfa
    b = 0.0                         # parametro calculado - beta
    fx = 0.0                        # valor da funcao distribuicao
    somaF = 0.0                     # somatorio dos valores da funcao de distribuicao
    result_ls = []                  # retorno: lista com parametros calculados [C1,a,b]
    
    N = len(faixa_ls)
    
    faixa_data = faixa_ls[0]
    x_min = faixa_data[1]
    
    faixa_data = faixa_ls[N - 1]
    x_max = faixa_data[2]
    
    dif = x_max - x_min
    
    # EXEC
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_01 = (faixa_meio - x_min) / dif
        faixa_prob = faixa_data[4]; 
        valor_exp = valor_exp + (faixa_01 * faixa_prob)
        valor_exp2 = valor_exp2 + (faixa_01 * faixa_01 * faixa_prob)

    valor_var = valor_exp2 - (valor_exp * valor_exp)

    valor_v = (valor_exp * (1.0 - valor_exp)) / valor_var

    # Distribuicao Beta 
    
    a = valor_exp * valor_v;
    b = (1 - valor_exp) * valor_v
        
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_01 = (faixa_meio - x_min) / dif
        fx = peam_func_pdf_beta(faixa_01, 1.0, a, b)
        somaF = somaF + fx
    
    C1 = 1.0 / somaF
    
    result_ls.append(C1)
    result_ls.append(a)
    result_ls.append(b)
    
    # RET
    return result_ls

# peam_mle_exp(): calcula o valor do parametro "lambda" em funcao dos valores das frequencias das faixas
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# retorno: lista com valores calculados dos parametros [C1,lambda]
def peam_mle_exp(faixa_ls):
    # VARS DEF
    #
    faixa_data = []                 # elemento da lista de faixas [Nfaixa,min,max,freq]
    faixa_meio = 0.0                # valor central da faixa
    faixa_prob = 0.0                # valor central da faixa
    valor_exp = 0.0                 # valor experado    
    val_lambda = 0.0                # parametro calculado - lambda
    C1 = 0.0                        # parametro calculado - C1
    fx = 0.0                        # valor da funcao distribuicao
    somaF = 0.0                     # somatorio dos valores da funcao de distribuicao
    result_ls = []                  # retorn: lista com parametros calculados [C1,lambda]
    
    # EXEC
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_prob = faixa_data[4]; 
        valor_exp = valor_exp + (faixa_meio * faixa_prob)

    # Distribuicao Exponencial 
    # E(X) = 1/lambda
    # VAR(X) = 1/(lambda * lambda)
    
    val_lambda = 1 / valor_exp;

    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        fx = peam_func_pdf_exp(faixa_meio, 1.0, val_lambda)
        somaF = somaF + fx
    
    C1 = 1.0 / somaF
    
    result_ls.append(C1)
    result_ls.append(val_lambda)
    
    # RET
    return result_ls

# peam_mle_norm(): calcula o valor do parametro "std" em funcao dos valores das frequencias das faixas
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# media - media dos valores
# retorno: lista contendo a media, o desvio padrao e a variancia calculados [media,std,var]
def peam_mle_norm(faixa_ls, media):
    # VARS DEF
    #
    faixa_data = []                 # elemento da lista de faixas [Nfaixa,min,max,freq]
    faixa_meio = 0.0                # valor central da faixa
    faixa_freq = 0.0                # valor da frequencia da faixa
    dist_media2 = 0.0               # quadrado da distancia ate o valor medio dos dados
    N = 0.0                         # numero de faixas
    val_exp = 0.0                   # valor calculado da media da distribuicao    
    val_var = 0.0                   # valor calculado da variancia da distribuicao
    val_std = 0.0                   # valor calculado do desvio padrao da distribuicao
    C1 = 0.0                        # valor calculado da constante de ajuste
    fx = 0.0                        # valor da funcao distribuicao
    somaF = 0.0                     # somatorio dos valores da funcao de distribuicao
    result_ls = []                  # retorno: parametros calculados - media e variancia
    
    # EXEC
    
    # calcula media e variancia da distribuicao
    #
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_prob = faixa_data[4] 
        dist_media2 = (faixa_meio - media) * (faixa_meio - media)
        val_exp = val_exp + (faixa_meio * faixa_prob)
        val_var = val_var + (dist_media2 * faixa_prob)

    val_std = M.sqrt(val_var)

    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        fx = peam_func_pdf_norm(faixa_meio, 1.0, val_exp, val_std)
        somaF = somaF + fx
    
    C1 = 1.0 / somaF
    
    result_ls.append(C1)
    result_ls.append(val_exp);
    result_ls.append(val_std);
    result_ls.append(val_var);
    
    # RET
    return result_ls

# peam_mle_lognormOLD(): calcula o valor do parametro "std" em funcao dos valores das frequencias das faixas
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# media - media dos valores
# retorno: lista contendo a media, o desvio padrao e a variancia calculados [media,std,var]
def peam_mle_lognormOLD(faixa_ls, media):
    # VARS DEF
    #
    faixa_data = []                 # elemento da lista de faixas [Nfaixa,min,max,freq]
    faixa_meio = 0.0                # valor central da faixa
    dist_media2 = 0.0               # quadrado da distancia ate o valor medio dos dados
    val_exp = 0.0                   # valor calculado da media da distribuicao    
    val_var = 0.0                   # valor calculado da variancia da distribuicao
    val_std = 0.0                   # valor calculado do desvio padrao da distribuicao
    C1 = 0.0                        # valor calculado da constante de ajuste
    fx = 0.0                        # valor da funcao distribuicao
    somaF = 0.0                     # somatorio dos valores da funcao de distribuicao
    result_ls = []                  # retorno: parametros calculados - media e variancia
    
    # EXEC
    
    # calcula media e variancia da distribuicao
    #
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_prob = faixa_data[4] 
        dist_media2 = (M.log(faixa_meio) - M.log(media)) * (M.log(faixa_meio) - M.log(media))
        val_exp = val_exp + (M.log(faixa_meio) * faixa_prob)
        val_var = val_var + (dist_media2 * faixa_prob)

    val_std = M.sqrt(val_var)

    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        fx = peam_func_pdf_lognorm(faixa_meio, 1.0, val_exp, val_std)
        somaF = somaF + fx
    
    C1 = 1.0 / somaF
    
    result_ls.append(C1)
    result_ls.append(val_exp);
    result_ls.append(val_std);
    result_ls.append(val_var);
    
    # RET
    return result_ls

# peam_mle_lognorm(): calcula o valor do parametro "std" em funcao dos valores das frequencias das faixas
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# media - media dos valores
# retorno: lista contendo a media, o desvio padrao e a variancia calculados [media,std,var]
def peam_mle_lognorm(faixa_ls, media):
    # VARS DEF
    #
    faixa_data = []                 # elemento da lista de faixas [Nfaixa,min,max,freq]
    faixa_meio = 0.0                # valor central da faixa
    dist_media2 = 0.0               # quadrado da distancia ate o valor medio dos dados
    val_exp = 0.0                   # valor calculado da media da distribuicao    
    val_var = 0.0                   # valor calculado da variancia da distribuicao
    val_std = 0.0                   # valor calculado do desvio padrao da distribuicao
    C1 = 0.0                        # valor calculado da constante de ajuste
    fx = 0.0                        # valor da funcao distribuicao
    somaF = 0.0                     # somatorio dos valores da funcao de distribuicao
    result_ls = []                  # retorno: parametros calculados - media e variancia
    
    # EXEC
    
    # calcula media e variancia da distribuicao
    #
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_prob = faixa_data[4] 
        val_exp = val_exp + (M.log(faixa_meio) * faixa_prob)

    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_prob = faixa_data[4] 
        dist_media2 = (M.log(faixa_meio) - val_exp) * (M.log(faixa_meio) - val_exp)
        val_var = val_var + (dist_media2 * faixa_prob)

    val_std = M.sqrt(val_var)

    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        fx = peam_func_pdf_lognorm(faixa_meio, 1.0, val_exp, val_std)
        somaF = somaF + fx
    
    C1 = 1.0 / somaF
    
    result_ls.append(C1)
    result_ls.append(val_exp);
    result_ls.append(val_std);
    result_ls.append(val_var);
    
    # RET
    return result_ls

# peam_mle_weibull_deriv(): calcula o valor da derivada da PDF da Weibull em funcao de "k"
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# k - valor de "k" considerado
# retorno: valor da derivada da PDF da Weibull
def peam_mle_weibull_deriv_part1(faixa_ls, k):
    # VARS DEF
    #
    sN1 = 0.0
    sD1 = 0.0
    s2 = 0.0
    val = 0.0
    n = float(len(faixa_ls))

    # EXEC
    
    # calcula "sN1" numerador da primeira parte da derivada
    # calcula "sD1" denominador da primeira parte da derivada
    # calcula "s2" segunda parte da derivada
    #
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        
        lnX = M.log(faixa_meio)
        
        xK = M.pow(faixa_meio, k)
        
        k_1 = 1.0 / k
        
        sN1 = sN1 + xK * lnX
        sD1 = sD1 + xK

        s2 = lnX

    if abs(sD1) > 1.0e-200:
        val = (sN1 / sD1) - k_1 - (s2 / n)
    else:
        s.exit(1)

    # RET
    return val

# peam_mle_weibull_deriv_part(): calcula o valor da derivada da PDF da Weibull em funcao de "k" por parte
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# k - valor de "k" considerado
# retorno: valor do valor das partes da derivada da PDF da Weibull [sN1,sD1,k_1,s2]
def peam_mle_weibull_deriv(faixa_ls, k):
    # VARS DEF
    #
    sN1 = 0.0
    sD1 = 0.0
    s2 = 0.0
    n = float(len(faixa_ls))
    result_ls = []

    # EXEC
    
    # calcula "sN1" numerador da primeira parte da derivada
    # calcula "sD1" denominador da primeira parte da derivada
    # calcula "s2" segunda parte da derivada
    #
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        
        lnX = M.log(faixa_meio)
        
        xK = M.pow(faixa_meio, k)
        
        k_1 = 1.0 / k
        
        sN1 = sN1 + xK * lnX
        sD1 = sD1 + xK

        s2 = lnX

    result_ls.append(sN1)
    result_ls.append(sD1)
    result_ls.append(- k_1)
    result_ls.append(- (s2 / n))

    # RET
    return result_ls

# peam_mle_weibull_k(): calcula o valor do parametro "k" variando o valor entre [-3.999999,9.999999]
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# retorno: lista contendo o valor estimado de "k" e do erro de aproximacao [k,err]
def peam_mle_weibull_k(faixa_ls):
    # VARS DEF
    #
    prec = 0.1
    min1 = -999999.999999
    max1 = -7.999999
    k = max1
    dx = 0.01
    val_ls = []
    result_ls = []

    # EXEC
    while k > min1:    
        val_ls = peam_mle_weibull_deriv(faixa_ls, k)

        sN1 = val_ls[0]
        sD1 = val_ls[1]
        k_1 = val_ls[2]
        s2 = val_ls[3]

        print("k=", k, "sN1=", sN1, "sD1=", sD1, "k_1=", k_1, "s2=", s2)
        
        if (abs(sN1) <= prec) and (abs(k_1) <= prec) and (abs(s2) <= prec):
            result_ls.append(k)
            result_ls.append(0.0)
            return k
        k = k - dx
        
    result_ls.append(k)

    # RET
    return result_ls

# peam_mle_weibull(): calcula o valor do parametro "lambda" em funcao de "k" e dos valores das frequencias das faixas
# faixa_ls - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
# k - valor de forma da curva
# retorno: lista contendo os valores [C1,lambda,k]
def peam_mle_weibull(faixa_ls, k):
    # VARS DEF
    #
    k_1 = 1.0 / k                   # valor inverso do parametro "k"
    faixa_data = []                 # elemento da lista de faixas [Nfaixa,min,max,freq]
    faixa_meio = 0.0                # valor central da faixa
    val_lambdaK = 0.0               # valor calculado do parametro "lambda" elevado a "k" da distribuicao    
    val_lambda = 0.0                # valor calculado do parametro "lambda" da distribuicao    
    C1 = 0.0                        # valor calculado da constante de ajuste
    fx = 0.0                        # valor da funcao distribuicao
    somaF = 0.0                     # somatorio dos valores da funcao de distribuicao
    result_ls = []                  # retorno: parametros calculados - media e variancia
    
    # EXEC
    
    # calcula media e variancia da distribuicao
    #
    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        faixa_prob = faixa_data[4]
        faixa_dist = M.pow(faixa_meio, k)
        val_lambdaK = val_lambdaK + (faixa_dist * faixa_prob)

    val_lambda = M.pow(val_lambdaK, k_1)

    for faixa_data in faixa_ls:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        fx = peam_func_pdf_weibull(faixa_meio, 1.0, val_lambda, k)
        somaF = somaF + fx
    
    C1 = 1.0 / somaF
    
    result_ls.append(C1)
    result_ls.append(val_lambda);
    result_ls.append(k);
    
    # RET
    return result_ls

# TESTES DE COMPARACAO DAS AMOSTRAS
#

# testOneWayANOVA
#

# testOneWayANOVA(): funcao que realiza o teste One-Way ANOVA para verificar se as medias das amostras sao iguais
# media1 - media da amostra 1
# std1 - desvio padrao da amostra 1
# media2 - media da amostra 2
# std2 - desvio padrao da amostra 2
# mediaC - media da amostra combinada
# stdC - desvio padrao da amostra combinada
# n - numero de grupos de populacao (parametro ignorado, default=2)
# m - numero de amostras em cada grupo (numero de dias =N)                 
# retorno: lista de valores contendo [w, MSb, MSw]
def peam_testOneWayANOVA(media1, std1, media2, std2, mediaC, stdC, n, m):
    # VARS DEF
    #
    n = 2
    #
    var1 = std1 * std1                  # variancia da amostra 1
    var2 = std2 * std2                  # variancia da amostra 2
    dist1_2 = (media1 - mediaC) * (media1 - mediaC)     # quadrado da distancia entre as medias da amostra 1 e do conjunto
    dist2_2 = (media2 - mediaC) * (media2 - mediaC)     # quadrado da distancia entre as medias da amostra 2 e do conjunto
    soma2 = dist1_2 + dist2_2           # soma dos quadrados das distancias das medias das amostras ao conjunto
    MSb = (m / (n - 1.0)) * soma2       # valor calculado da variancia entre grupos
    MSw = (var1 + var2) / n             # valor calculado da media das variancias dos grupos
    w = MSb / MSw                       # w-statistic
    result_ls = []                      # retorno: lista de resultado do teste [w, MSb, MSw]

    result_ls.append(w)
    result_ls.append(MSb)
    result_ls.append(MSw)
    
    # RET
    return result_ls

# testChiSquare
#

# testChiSquare(): funcao que realiza o teste Chi-Square para verificar se as medias da amostra aproxima uma distribuicao
# tfaixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1,f1,E1],[faixa2,min2,max2,freq2,p2,f2,E2], ... ,[faixaN,minN,maxN,freqN,pN,fN,EN] ]
# retorno: lista de valores contendo [G, X2]
def peam_testChiSquare(tfaixa_ls1):
    # VARS DEF
    #
    prec = 0.0000001
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq]
    faixa_meio1 = 0.0               # valor central da faixa
    result_ls = []                  # retorno: lista de resultado do teste [G, X2]
    somaG = 0.0
    somaX2 = 0.0
    dG = 0.0
    dX2 = 0.0
    Oi = 0.0
    Freq_i = 0.0
    Pi = 0.0
    fi = 0.0
    Ei = 0.0
    FreqN = 0.0

    for faixa_data1 in tfaixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0

        Oi = float(faixa_data1[0])
        Freq_i = faixa_data1[3]
        Pi = faixa_data1[4]
        fi = faixa_data1[5]
        Ei = faixa_data1[6]
        FreqN = faixa_data1[7]
        
        if Oi > prec and Ei > prec:
            dG = Oi * M.log(Oi / Ei)
        else:
            dG = 0.0
        
        if Ei > prec:
            dX2 = ((Oi - Ei) * (Oi - Ei)) / Ei
        else:
            dX2 = 0.0

        somaG = somaG + dG
        somaX2 = somaX2 + dX2

    val_G = 2 * somaG
    val_X2 = somaX2

    result_ls.append(val_G)
    result_ls.append(val_X2)
    
    # RET
    return result_ls

# testKS - Kolmogorov-Smirnov
#

# peam_testKS_exp(): funcao que realiza o teste Kolmogorov-Smirnov para verificar se as medias da amostra aproxima uma distribuicao
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# val_lambda_reg - valor de lambda regularizado para o intervalo [0,1] (lambda = 1.0 / media * (max_interval - min_interval))
# retorno: retorna o valor Dn calculado (distancia maxima entre a distribuicao empirica e a distribuicoes estimada)
def peam_testKS_exp(faixa_ls1, val_lambda_reg):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq]
    tmp_val = 0.0                   # acumulador para a probabilidade
    num_faixa1 = 0.0                # numero da faixa atual
    num_faixa_reg1 = 0.0            # numero da faixa proporcional ao intervalo [0,1]
    Fx = 0.0                        # valor da funcao CDF no ponto
    Dn = 0.0                        # retorno: valor calculado Dn

    Nfaixas = float(len(faixa_ls1)) # numero de faixas na lista

    num_faixa1 = 0.0
    for faixa_data1 in faixa_ls1:
        num_faixa_reg1 = num_faixa1 / Nfaixas

        Fx = peam_func_cdf_exp(num_faixa_reg1, 1.0, val_lambda_reg)

        tmp_val = tmp_val + faixa_data1[4]      # acumulador da probabilidade

        tmp_Dn = abs(tmp_val - Fx)        
        if tmp_Dn > Dn:
            Dn = tmp_Dn
    
    # RET
    return Dn

# peam_testKS_norm(): funcao que realiza o teste Kolmogorov-Smirnov para verificar se as medias da amostra aproxima uma distribuicao
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# media_reg - valor medio regularizado [0,1]
# std_reg - desvio padrao regularizado [0,1]
# retorno: retorna o valor Dn calculado (distancia maxima entre a distribuicao empirica e a distribuicoes estimada)
def peam_testKS_norm(faixa_ls1, media_reg, std_reg):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq]
    tmp_val = 0.0                   # acumulador para a probabilidade
    num_faixa1 = 0.0                # numero da faixa atual
    num_faixa_reg1 = 0.0            # numero da faixa proporcional ao intervalo [0,1]
    Fx = 0.0                        # valor da funcao CDF no ponto
    Dn = 0.0                        # retorno: valor calculado Dn

    Nfaixas = float(len(faixa_ls1)) # numero de faixas na lista

    num_faixa1 = 0.0
    for faixa_data1 in faixa_ls1:
        num_faixa_reg1 = num_faixa1 / Nfaixas

        Fx = peam_func_cdf_norm(num_faixa_reg1, 1.0, media_reg, std_reg)

        tmp_val = tmp_val + faixa_data1[4]      # acumulador da probabilidade

        tmp_Dn = abs(tmp_val - Fx)        
        if tmp_Dn > Dn:
            Dn = tmp_Dn
    
    # RET
    return Dn

# peam_testKS_lognorm(): funcao que realiza o teste Kolmogorov-Smirnov para verificar se as medias da amostra aproxima uma distribuicao
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# media_reg - valor medio regularizado [0,1]
# std_reg - desvio padrao regularizado [0,1]
# retorno: retorna o valor Dn calculado (distancia maxima entre a distribuicao empirica e a distribuicoes estimada)
def peam_testKS_lognorm(faixa_ls1, media_reg, std_reg):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq]
    tmp_val = 0.0                   # acumulador para a probabilidade
    num_faixa1 = 0.0                # numero da faixa atual
    num_faixa_reg1 = 0.0            # numero da faixa proporcional ao intervalo [0,1]
    Fx = 0.0                        # valor da funcao CDF no ponto
    Dn = 0.0                        # retorno: valor calculado Dn

    Nfaixas = float(len(faixa_ls1)) # numero de faixas na lista

    num_faixa1 = 0.0
    for faixa_data1 in faixa_ls1:
        num_faixa_reg1 = num_faixa1 / Nfaixas

        Fx = peam_func_cdf_lognorm(num_faixa_reg1, 1.0, media_reg, std_reg)

        tmp_val = tmp_val + faixa_data1[4]      # acumulador da probabilidade

        tmp_Dn = abs(tmp_val - Fx)        
        if tmp_Dn > Dn:
            Dn = tmp_Dn
    
    # RET
    return Dn

# peam_testKS_weibull(): funcao que realiza o teste Kolmogorov-Smirnov para verificar se as medias da amostra aproxima uma distribuicao
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# val_lambda_reg - valor de lambda regularizado para o intervalo [0,1] (lambda = 1.0 / media * (max_interval - min_interval))
# k - valor de "k"
# retorno: retorna o valor Dn calculado (distancia maxima entre a distribuicao empirica e a distribuicoes estimada)
def peam_testKS_weibull(faixa_ls1, val_lambda_reg, k):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq]
    tmp_val = 0.0                   # acumulador para a probabilidade
    num_faixa1 = 0.0                # numero da faixa atual
    num_faixa_reg1 = 0.0            # numero da faixa proporcional ao intervalo [0,1]
    Fx = 0.0                        # valor da funcao CDF no ponto
    Dn = 0.0                        # retorno: valor calculado Dn

    Nfaixas = float(len(faixa_ls1)) # numero de faixas na lista

    num_faixa1 = 0.0
    for faixa_data1 in faixa_ls1:
        num_faixa_reg1 = num_faixa1 / Nfaixas

        Fx = peam_func_cdf_weibull(num_faixa_reg1, 1.0, val_lambda_reg, k)

        tmp_val = tmp_val + faixa_data1[4]      # acumulador da probabilidade

        tmp_Dn = abs(tmp_val - Fx)        
        if tmp_Dn > Dn:
            Dn = tmp_Dn
    
    # RET
    return Dn

# CRIACAO DE GRAFICOS
#

# peam_show_boxplot(): apresenta histograma referente as frequencias das faixas
# titulo - titulo do grafico
# ds_ls11 - lista de elementos [ [key1,val1],[key2,val2], ... ,[keyN,valN] ]
# workset_ls1 - lista de elementos [ [key1,val1],[key2,val2], ... ,[keyN,valN] ]
def peam_show_boxplot(titulo, ds_ls1, workset_ls1):
    # VARS DEF
    #
    ds_data1 = []
    gkey_ls1 = []                   # lista de chaves do dataset
    gval_ls1 = []                   # lista de valores do dataset
    gkey_ls2 = []                   # lista de chaves do workset
    gval_ls2 = []                   # lista de valores do workset
        
    # EXEC
    for ds_data1 in ds_ls1:
        gkey_ls1.append(ds_data1[0])
        gval_ls1.append(ds_data1[1])

    for ds_data1 in workset_ls1:
        gkey_ls2.append(ds_data1[0])
        gval_ls2.append(ds_data1[1])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    str_title0 = "DATASET"
    axs[0].set_title(str_title0)
    axs[0].boxplot(gval_ls1, notch=True)

    str_title1 = "WORKSET"
    axs[1].set_title(str_title1)
    axs[1].boxplot(gval_ls2, notch=True)
    
    plt.show()

# peam_show_hist(): apresenta histograma referente as frequencias das faixas
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1],[faixa2,min2,max2,freq2], ... ,[faixaN,minN,maxN,freqN] ]
def peam_show_hist(titulo, faixa_ls1):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq]
    faixa_meio1 = 0.0               # valor central da faixa
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gfreq_ls1 = []                  # lista de frequencias
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gfreq_ls1.append(faixa_data1[3])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gfreq_ls1, width=1.0)
    axs[1].plot(gval_ls1, gfreq_ls1)
    
    plt.show()

# peam_show_erf(): apresenta o grafico da funcao de erro de Gauss
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
def peam_show_erf(titulo, faixa_ls1):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gF_ls1 = []                     # lista de valores da funcao
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_erf(faixa_meio1)        
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gF_ls1.append(fx)
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gF_ls1, width=1.0)
    axs[1].plot(gval_ls1, gF_ls1)
    
    plt.show()

# peam_show_fpm_geom(): apresenta o grafico da funcao probabilidade de massa da Distribuicao Geometrica
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# pS - probabilida de ocorrencia de sucesso em uma distribuicao geometrica
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_fpm_geom(titulo, faixa_ls1, C1, pS, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_fpm_geom(faixa_meio1, C1, pS)        
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(faixa_data1[4])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gP_ls1, width=1.0)
    if( bSobreposicao ):
        axs[0].bar(gfaixa_ls1, gP_ls2, width=1.0)
        

    axs[1].plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs[1].plot(gval_ls1, gP_ls2)
    
    plt.show()

# peam_show_fpa_geom(): apresenta o grafico da funcao probabilidade acumulada da Distribuicao Geometrica
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# pS - probabilida de ocorrencia de sucesso em uma distribuicao geometrica
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_fpa_geom(titulo, faixa_ls1, C1, pS, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    tmp_val1 = 0.0                  # acumula valores da distribuicao
    tmp_val2 = 0.0                  # acumula valores da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_fpm_geom(faixa_meio1, C1, pS)        
        
        tmp_val1 = tmp_val1 + fx     
        tmp_val2 = tmp_val2 + faixa_data1[4]     
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(tmp_val1)
        gP_ls2.append(tmp_val2)
        
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs.plot(gval_ls1, gP_ls2)
    
    plt.show()

# peam_show_qqplot_geom(): apresenta o grafico QQ-PLOT da funcao de Distribuicao Geometrica
# titulo - titulo da faixa de elementos
# ds1 - lista de elementos do dataset [ [key1,val1],[key2,val2], ... ,[keyN,valN] ]
# pS - probabilida de ocorrencia de sucesso em uma distribuicao geometrica
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_qqplot_geom(titulo, ds1, C1, pS, bSobreposicao, min1, max1, N1, div, scale):
    # VARS DEF
    #
    N_faixas = M.floor(N1 / div) - 1
    N_min = (N_faixas / 4.0) - 1
    N_max = N1 - (N_faixas / 4.0) + 1
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx_i1 = 0.0                     # valor amostral da funcao
    fx_f1 = 0.0                     # valor amostral da funcao
    fx_i2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    fx_f2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    gval_ls1 = []                   # lista de valores medios da faixa
    gR_ls1 = []                     # razao entre o valor calculado e o medido
    gR_ls2 = []                     # linha a 45g

    faixa_ls1 = peam_lista_freq(ds1, N_faixas, min1, max1, N1)        

    # EXEC
    n_i = 0.0
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0

        fx_f1 = faixa_data1[4]
        fx_f2 = peam_func_fpm_geom(faixa_meio1, C1, pS)        

        if ((n_i > N_min) and (n_i < N_max)):        
            dif1 = fx_f1 - fx_i1
            dif2 = fx_f2 - fx_i2

            slope = (dif1 / dif2) * scale

            x = n_i + 0.5            
            y1 = x - slope * x
            y2 = x

            #print("n_i:", n_i, "x:", x, "y1:", y1, "y2:", y2)

            gval_ls1.append(x)
            gR_ls1.append(y1)
            gR_ls2.append(y2)

        fx_i1 = fx_f1
        fx_i2 = fx_f2

        n_i = n_i + 1.0

                
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.scatter(gval_ls1, gR_ls1)
    axs.plot(gval_ls1, gR_ls2)
        
    plt.show()

# peam_show_pdf_beta(): apresenta o grafico da funcao densidade de distribuicao da Distribuicao Beta
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# C1 - constante de ajuste
# a - parametro "a" da distribuicao 
# b - parametro "b" da distribuicao
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_pdf_beta(titulo, faixa_ls1, C1, a, b, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    faixa_01 = 0.0                  # valor central da faixa regularizada [0,1]
    fx = 0.0                        # valor calculado para a funcao Beta
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
         
    N = len(faixa_ls1)
    
    faixa_data1 = faixa_ls1[0]
    x_min = faixa_data1[1]
    
    faixa_data1 = faixa_ls1[N - 1]
    x_max = faixa_data1[2]   
    
    dif = x_max - x_min
    
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        faixa_01 = (faixa_meio1 - x_min) / dif
        
        fx = peam_func_pdf_beta(faixa_01, C1, a, b)   
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_01)
        gP_ls1.append(fx)
        gP_ls2.append(faixa_data1[4])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gP_ls1, width=1.0)
    if( bSobreposicao ):
        axs[0].bar(gfaixa_ls1, gP_ls2, width=1.0)
        

    axs[1].plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs[1].plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_cdf_beta(): apresenta o grafico da funcao de distribuicao acumulada da Distribuicao Beta
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# C1 - constante de ajuste
# a - parametro "a" da distribuicao 
# b - parametro "b" da distribuicao
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_cdf_beta(titulo, faixa_ls1, C1, a, b, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    faixa_01 = 0.0                  # valor central da faixa regularizada [0,1]
    tmp_val = 0.0                   # acumula o valor da distribuicao
    fx = 0.0                        # valor calculado para a funcao Beta
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
         
    N = len(faixa_ls1)
    
    faixa_data1 = faixa_ls1[0]
    x_min = faixa_data1[1]
    
    faixa_data1 = faixa_ls1[N - 1]
    x_max = faixa_data1[2]   
    
    dif = x_max - x_min
    
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        faixa_01 = (faixa_meio1 - x_min) / dif
        
        fx = peam_func_cdf_beta(faixa_01, C1, a, b)   
        
        tmp_val = tmp_val + faixa_data1[4]
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_01)
        gP_ls1.append(fx)
        gP_ls2.append(tmp_val)
        
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs.plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_pdf_beta_01(): apresenta o grafico da funcao densidade de distribuicao da Distribuicao Beta aplicado na F
# d1 - parametro "d1" da distribuicao F 
# d2 - parametro "d2" da distribuicao F
def peam_show_pdf_beta_01(d1, d2):
    # VARS DEF
    #
    min_x = 0.01                    # valor minimo da variavel aleatoria
    max_x = 1.00                    # valor maximo da variavel aleatoria
    dx = 0.01                       # valor do incremento da variavel aleatoria
    x = min_x                       # valor da variavel aleatoria
    y = 0.0                         # valor calculado para a funcao F    
    a = d1 / 2.0                    # valor do parametro "alfa" da distribuicao Beta
    b = d2 / 2.0                    # valor do parametro "beta" da distribuicao Beta
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao

    while x <= max_x:
        y = peam_func_pdf_beta(x, 1.0, a, b)   
        gval_ls1.append(x)
        gP_ls1.append(y)
        x = x + dx
        
    str_titulo = "Beta(x,alfa=" + str(d1) + ",beta=" + str(d2) + ") no intervalo [0,1]"
        
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(str_titulo)
    
    axs.plot(gval_ls1, gP_ls1)
        
    plt.show()

# peam_show_pdf_F(): apresenta o grafico da funcao densidade de distribuicao da Distribuicao F (Fisher-Snedecor)
# d1 - parametro "d1" da distribuicao 
# d2 - parametro "d2" da distribuicao
# scale - escala para os valores de "x"
# x_conf - posicao da marcacao de "x"
# x_pval - posicao da marcacao do "p-value" 
def peam_show_pdf_F(d1, d2, scale, x_conf, x_pval):
    # VARS DEF
    #
    min_x = 0.01                     # valor minimo da variavel aleatoria
    max_x = 99.99                    # valor maximo da variavel aleatoria
    dx = 0.01                        # valor do incremento da variavel aleatoria
    x = min_x                        # valor da variavel aleatoria
    x_scale = x * scale              # valor da variavel aleatoria
    y = 0.0                          # valor calculado para a funcao F    
    gval_ls1 = []                    # lista de valores medios da faixa
    gP_ls1 = []                      # lista de valores da funcao distribuicao

    while x <= max_x:
        x_scale = x * scale
        y = peam_func_pdf_F(x_scale, d1, d2)   
        gval_ls1.append(x_scale)
        gP_ls1.append(y)
        x = x + dx
        
    str_titulo = "pdf-F(x,d1=" + str(d1) + ",d2=" + str(d2) + ") no intervalo [0,1]"

    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(str_titulo)
    
    axs.plot(gval_ls1, gP_ls1)

    gloc = []
    gloc.append(x_conf)
    gloc.append(x_pval)

    glbl = []
    glbl.append("x: confidence")
    glbl.append("x: p-value")    

    axs.set_xticks(gloc)
    axs.set_xticklabels(glbl)
        
    plt.show()

# peam_show_cdf_F(): apresenta o grafico da funcao distribuicao acumulada da Distribuicao F (Fisher-Snedecor)
# d1 - parametro "d1" da distribuicao 
# d2 - parametro "d2" da distribuicao
# scale - escala para os valores de "x"
# x_conf - posicao da marcacao de "x"
# x_pval - posicao da marcacao do "p-value" 
def peam_show_cdf_F(d1, d2, scale, x_conf, x_pval):
    # VARS DEF
    #
    min_x = 0.01                     # valor minimo da variavel aleatoria
    max_x = 99.99                    # valor maximo da variavel aleatoria
    dx = 0.01                        # valor do incremento da variavel aleatoria
    x = min_x                        # valor da variavel aleatoria
    x_scale = x * scale              # valor da variavel aleatoria
    y = 0.0                          # valor calculado para a funcao F    
    gval_ls1 = []                    # lista de valores medios da faixa
    gP_ls1 = []                      # lista de valores da funcao distribuicao

    while x <= max_x:
        x_scale = x * scale
        y = peam_func_cdf_F(x_scale, d1, d2)   
        gval_ls1.append(x_scale)
        gP_ls1.append(y)
        x = x + dx
        
    str_titulo = "cdf-F(x,d1=" + str(d1) + ",d2=" + str(d2) + ") no intervalo [0,1]"

    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(str_titulo)
    
    axs.plot(gval_ls1, gP_ls1)

    gloc = []
    gloc.append(x_conf)
    gloc.append(x_pval)

    glbl = []
    glbl.append("x: confidence")
    glbl.append("x: p-value")    

    axs.set_xticks(gloc)
    axs.set_xticklabels(glbl)
        
    plt.show()

# peam_show_pdf_exp(): apresenta o grafico da funcao densidade de distribuicao da Distribuicao Exponencial
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# C1 - contante de ajuste dos dados da exponencial
# val_lambda - media de ocorrencia dos eventos
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_pdf_exp(titulo, faixa_ls1, C1, val_lambda, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_pdf_exp(faixa_meio1, C1, val_lambda)        
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(faixa_data1[4])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gP_ls1, width=1.0)
    if( bSobreposicao ):
        axs[0].bar(gfaixa_ls1, gP_ls2, width=1.0)
        

    axs[1].plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs[1].plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_cdf_exp(): apresenta o grafico da funcao de distribuicao acumulada da Distribuicao Exponencial
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# C1 - contante de ajuste dos dados da exponencial
# val_lambda - media de ocorrencia dos eventos
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_cdf_exp(titulo, faixa_ls1, C1, val_lambda, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    tmp_val = 0.0                   # valor acumulado dos elementos da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_cdf_exp(faixa_meio1, C1, val_lambda)        
        
        tmp_val = tmp_val + faixa_data1[4]

        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(tmp_val)
        
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs.plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_qqplot_exp(): apresenta o grafico QQ-PLOT da funcao densidade de distribuicao da Distribuicao Exponencial
# titulo - titulo da faixa de elementos
# ds1 - lista de elementos do dataset [ [key1,value1], ... ,[keyN,valueN] ]
# C1 - contante de ajuste dos dados da exponencial
# val_lambda - media de ocorrencia dos eventos
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_qqplot_exp(titulo, ds1, C1, val_lambda, bSobreposicao, min1, max1, N1, div, scale):
    # VARS DEF
    #
    N_faixas = M.floor(N1 / div) - 1
    N_min = (N_faixas / 4.0) - 1
    N_max = N1 - (N_faixas / 4.0) + 1
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx_i1 = 0.0                     # valor amostral da funcao
    fx_f1 = 0.0                     # valor amostral da funcao
    fx_i2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    fx_f2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    gval_ls1 = []                   # lista de valores medios da faixa
    gR_ls1 = []                     # razao entre o valor calculado e o medido
    gR_ls2 = []                     # linha a 45g

    faixa_ls1 = peam_lista_freq(ds1, N_faixas, min1, max1, N1)        

    # EXEC
    n_i = 0.0
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0

        fx_f1 = faixa_data1[4]
        fx_f2 = peam_func_pdf_exp(faixa_meio1, C1, val_lambda)        

        if ((n_i > N_min) and (n_i < N_max)):        
            dif1 = fx_f1 - fx_i1
            dif2 = fx_f2 - fx_i2

            slope = (dif1 / dif2) * scale

            x = n_i + 0.5            
            y1 = x - slope * x
            y2 = x

            #print("n_i:", n_i, "x:", x, "y1:", y1, "y2:", y2)

            gval_ls1.append(x)
            gR_ls1.append(y1)
            gR_ls2.append(y2)

        fx_i1 = fx_f1
        fx_i2 = fx_f2

        n_i = n_i + 1.0

                
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.scatter(gval_ls1, gR_ls1)
    axs.plot(gval_ls1, gR_ls2)
        
    plt.show()

# peam_show_pdf_norm(): apresenta o grafico da funcao densidade de distribuicao da Distribuicao Normal
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# media - media da distribuicao
# std - desvio padrao da distribuicao
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_pdf_norm(titulo, faixa_ls1, C1, media, std, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_pdf_norm(faixa_meio1, C1, media, std)        
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(faixa_data1[4])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gP_ls1, width=1.0)
    if( bSobreposicao ):
        axs[0].bar(gfaixa_ls1, gP_ls2, width=1.0)
        
    axs[1].plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs[1].plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_cdf_norm(): apresenta o grafico da funcao de distribuicao acumulada da Distribuicao Normal
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# media - media da distribuicao
# std - desvio padrao da distribuicao
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_cdf_norm(titulo, faixa_ls1, C1, media, std, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    tmp_val = 0.0                   # valor acumulado dos elementos da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_cdf_norm(faixa_meio1, C1, media, std)        

        tmp_val = tmp_val + faixa_data1[4]
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(tmp_val)
        
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs.plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_qqplot_norm(): apresenta o grafico QQ-PLOT da funcao densidade de distribuicao da Distribuicao Normal
# titulo - titulo da faixa de elementos
# ds1 - lista de elementos do dataset [ [key1,value1], ... ,[keyN,valueN] ]
# C1 - contante de ajuste dos dados da exponencial
# media - media da distribuicao
# std - desvio padrao da distribuicao
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_qqplot_norm(titulo, ds1, C1, media, std, bSobreposicao, min1, max1, N1, div, scale):
    # VARS DEF
    #
    N_faixas = M.floor(N1 / div) - 1
    N_min = (N_faixas / 4.0) - 1
    N_max = N1 - (N_faixas / 4.0) + 1
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx_i1 = 0.0                     # valor amostral da funcao
    fx_f1 = 0.0                     # valor amostral da funcao
    fx_i2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    fx_f2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    gval_ls1 = []                   # lista de valores medios da faixa
    gR_ls1 = []                     # razao entre o valor calculado e o medido
    gR_ls2 = []                     # linha a 45g

    faixa_ls1 = peam_lista_freq(ds1, N_faixas, min1, max1, N1)        

    # EXEC
    n_i = 0.0
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0

        fx_f1 = faixa_data1[4]
        fx_f2 = peam_func_pdf_norm(faixa_meio1, C1, media, std)        

        if ((n_i > N_min) and (n_i < N_max)):        
            dif1 = fx_f1 - fx_i1
            dif2 = fx_f2 - fx_i2

            slope = (dif1 / dif2) * scale

            x = n_i + 0.5            
            y1 = x - slope * x
            y2 = x

            #print("n_i:", n_i, "x:", x, "y1:", y1, "y2:", y2)

            gval_ls1.append(x)
            gR_ls1.append(y1)
            gR_ls2.append(y2)

        fx_i1 = fx_f1
        fx_i2 = fx_f2

        n_i = n_i + 1.0

                
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.scatter(gval_ls1, gR_ls1)
    axs.plot(gval_ls1, gR_ls2)
        
    plt.show()

# peam_show_pdf_lognorm(): apresenta o grafico da funcao densidade de distribuicao da Distribuicao Log-normal
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# media - media da distribuicao
# std - desvio padrao da distribuicao
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_pdf_lognorm(titulo, faixa_ls1, C1, media, std, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx = 0.0                        # valor calculado para a funcao
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_pdf_lognorm(faixa_meio1, C1, media, std)        
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(faixa_data1[4])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gP_ls1, width=1.0)
    if( bSobreposicao ):
        axs[0].bar(gfaixa_ls1, gP_ls2, width=1.0)
        
    axs[1].plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs[1].plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_cdf_lognorm(): apresenta o grafico da funcao de distribuicao acumulada da Distribuicao Log-normal
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# media - media da distribuicao
# std - desvio padrao da distribuicao
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_cdf_lognorm(titulo, faixa_ls1, C1, media, std, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    tmp_val = 0.0                   # valor acumulado dos elementos da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_cdf_lognorm(faixa_meio1, C1, media, std)        

        tmp_val = tmp_val + faixa_data1[4]
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(tmp_val)
        
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs.plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_qqplot_lognorm(): apresenta o grafico QQ-PLOT da funcao densidade de distribuicao da Distribuicao Log-normal
# titulo - titulo da faixa de elementos
# ds1 - lista de elementos do dataset [ [key1,value1], ... ,[keyN,valueN] ]
# C1 - contante de ajuste dos dados da exponencial
# media - media da distribuicao
# std - desvio padrao da distribuicao
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_qqplot_lognorm(titulo, ds1, C1, media, std, bSobreposicao, min1, max1, N1, div, scale):
    # VARS DEF
    #
    N_faixas = M.floor(N1 / div) - 1
    N_min = (N_faixas / 4.0) - 1
    N_max = N1 - (N_faixas / 4.0) + 1
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx_i1 = 0.0                     # valor amostral da funcao
    fx_f1 = 0.0                     # valor amostral da funcao
    fx_i2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    fx_f2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    gval_ls1 = []                   # lista de valores medios da faixa
    gR_ls1 = []                     # razao entre o valor calculado e o medido
    gR_ls2 = []                     # linha a 45g

    faixa_ls1 = peam_lista_freq(ds1, N_faixas, min1, max1, N1)        

    # EXEC
    n_i = 0.0
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0

        fx_f1 = faixa_data1[4]
        fx_f2 = peam_func_pdf_lognorm(faixa_meio1, C1, media, std)        

        if ((n_i > N_min) and (n_i < N_max)):        
            dif1 = fx_f1 - fx_i1
            dif2 = fx_f2 - fx_i2

            slope = (dif1 / dif2) * scale

            x = n_i + 0.5            
            y1 = x - slope * x
            y2 = x

            #print("n_i:", n_i, "x:", x, "y1:", y1, "y2:", y2)

            gval_ls1.append(x)
            gR_ls1.append(y1)
            gR_ls2.append(y2)

        fx_i1 = fx_f1
        fx_i2 = fx_f2

        n_i = n_i + 1.0

                
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.scatter(gval_ls1, gR_ls1)
    axs.plot(gval_ls1, gR_ls2)
        
    plt.show()

# peam_show_pdf_weibull(): apresenta o grafico da funcao densidade de distribuicao da Distribuicao Weibull
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# C1 - contante de ajuste dos dados da exponencial
# val_lambda - valor do parametro "lambda" da distribuicao
# k - valor do parametro "k" que define a forma da curva
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_pdf_weibull(titulo, faixa_ls1, C1, val_lambda, k, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_pdf_weibull(faixa_meio1, C1, val_lambda, k)        
        
        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(faixa_data1[4])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gP_ls1, width=1.0)
    if( bSobreposicao ):
        axs[0].bar(gfaixa_ls1, gP_ls2, width=1.0)
        
    axs[1].plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs[1].plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_cdf_weibull(): apresenta o grafico da funcao de distribuicao acumulada da Distribuicao Weibull
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# C1 - contante de ajuste dos dados da exponencial
# val_lambda - valor do parametro "lambda" da distribuicao
# k - valor do parametro "k" que define a forma da curva
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_cdf_weibull(titulo, faixa_ls1, C1, val_lambda, k, bSobreposicao):
    # VARS DEF
    #
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    tmp_val = 0.0                   # valor acumulado dos elementos da faixa
    fx = 0.0                        # valor calculado para a funcao geometrica dado um valor de "p"
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0
        
        fx = peam_func_cdf_weibull(faixa_meio1, C1, val_lambda, k)        
        
        tmp_val = tmp_val + faixa_data1[4]

        gfaixa_ls1.append(faixa_data1[0])
        gval_ls1.append(faixa_meio1)
        gP_ls1.append(fx)
        gP_ls2.append(tmp_val)
        
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.plot(gval_ls1, gP_ls1)
    if( bSobreposicao ):
        axs.plot(gval_ls1, gP_ls2)
        
    plt.show()

# peam_show_qqplot_weibull(): apresenta o grafico QQ-PLOT da funcao densidade de distribuicao da Distribuicao Weibull
# titulo - titulo da faixa de elementos
# ds1 - lista de elementos do dataset [ [key1,value1], ... ,[keyN,valueN] ]
# C1 - contante de ajuste dos dados da exponencial
# val_lambda - valor do parametro "lambda" da distribuicao
# k - valor do parametro "k" que define a forma da curva
# bSobreposicao - flag de sobreposicao das medicoes
def peam_show_qqplot_weibull(titulo, ds1, C1, val_lambda, k, bSobreposicao, min1, max1, N1, div, scale):
    # VARS DEF
    #
    N_faixas = M.floor(N1 / div) - 1
    N_min = (N_faixas / 4.0) - 1
    N_max = N1 - (N_faixas / 4.0) + 1
    faixa_data1 = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio1 = 0.0               # valor central da faixa
    fx_i1 = 0.0                     # valor amostral da funcao
    fx_f1 = 0.0                     # valor amostral da funcao
    fx_i2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    fx_f2 = 0.0                     # valor calculado para a funcao geometrica dado um valor de "p"
    gval_ls1 = []                   # lista de valores medios da faixa
    gR_ls1 = []                     # razao entre o valor calculado e o medido
    gR_ls2 = []                     # linha a 45g

    faixa_ls1 = peam_lista_freq(ds1, N_faixas, min1, max1, N1)        

    # EXEC
    n_i = 0.0
    for faixa_data1 in faixa_ls1:
        faixa_meio1 = faixa_data1[1] + (faixa_data1[2] - faixa_data1[1]) / 2.0

        fx_f1 = faixa_data1[4]
        fx_f2 = peam_func_pdf_weibull(faixa_meio1, C1, val_lambda, k)        

        if ((n_i > N_min) and (n_i < N_max)):        
            dif1 = fx_f1 - fx_i1
            dif2 = fx_f2 - fx_i2
            
            slope = (dif1 / dif2) * scale

            x = n_i + 0.5            
            y1 = x - slope * x
            y2 = x

            #print("n_i:", n_i, "x:", x, "y1:", y1, "y2:", y2)

            gval_ls1.append(x)
            gR_ls1.append(y1)
            gR_ls2.append(y2)

        fx_i1 = fx_f1
        fx_i2 = fx_f2

        n_i = n_i + 1.0

                
    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs.scatter(gval_ls1, gR_ls1)
    axs.plot(gval_ls1, gR_ls2)
        
    plt.show()

# peam_show_hist2(): apresenta o grafico do histograma de duas funcoes sobrepostas
# titulo - titulo da faixa de elementos
# faixa_ls1 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
# faixa_ls2 - lista de elementos por faixa [ [faixa1,min1,max1,freq1,p1],[faixa2,min2,max2,freq2,p2], ... ,[faixaN,minN,maxN,freqN,pN] ]
def peam_show_hist2(titulo, faixa_ls1, faixa_ls2):
    # VARS DEF
    #
    faixa_data = []                # elemento da lista de faixas [Nfaixa,min,max,freq,pN]
    faixa_meio = 0.0               # valor central da faixa
    gfaixa_ls1 = []                 # lista de faixas
    gval_ls1 = []                   # lista de valores medios da faixa
    gP_ls1 = []                     # lista de valores da funcao distribuicao
    gfaixa_ls2 = []                 # lista de faixas
    gval_ls2 = []                   # lista de valores medios da faixa
    gP_ls2 = []                     # lista de probabilidades para cada faixa
        
    # EXEC
    for faixa_data in faixa_ls1:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        
        gfaixa_ls1.append(faixa_data[0])
        gval_ls1.append(faixa_meio)
        gP_ls1.append(faixa_data[4])
        
    for faixa_data in faixa_ls2:
        faixa_meio = faixa_data[1] + (faixa_data[2] - faixa_data[1]) / 2.0
        
        gfaixa_ls2.append(faixa_data[0])
        gval_ls2.append(faixa_meio)
        gP_ls2.append(faixa_data[4])
        
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    fig.suptitle(titulo)
    
    axs[0].bar(gfaixa_ls1, gP_ls1, width=1.0)
    axs[0].bar(gfaixa_ls2, gP_ls2, width=1.0)
        
    axs[1].plot(gval_ls1, gP_ls1)
    axs[1].plot(gval_ls2, gP_ls2)
        
    plt.show()

# peam_show_pdf_chi2(): apresenta o grafico da funcao densidade de distribuicao da Distribuicao Chi-Square
# df - parametro "df" da distribuicao 
# scale - escala para os valores de "x"
# x_conf - posicao da marcacao de "x"
# x_pval - posicao da marcacao do "p-value" 
def peam_show_pdf_chi2(df, scale, x_conf, x_pval):
    # VARS DEF
    #
    min_x = 0.01                     # valor minimo da variavel aleatoria
    max_x = 99.99                    # valor maximo da variavel aleatoria
    dx = 0.01                        # valor do incremento da variavel aleatoria
    x = min_x                        # valor da variavel aleatoria
    x_scale = x * scale              # valor da variavel aleatoria
    y = 0.0                          # valor calculado para a funcao F    
    gval_ls1 = []                    # lista de valores medios da faixa
    gP_ls1 = []                      # lista de valores da funcao distribuicao

    while x <= max_x:
        x_scale = x * scale
        y = peam_func_pdf_chi2(x_scale, df)   
        gval_ls1.append(x_scale)
        gP_ls1.append(y)
        x = x + dx
        
    str_titulo = "pdf-CHI2(x,df=" + str(df) + ") no intervalo [0,1]"

    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(str_titulo)
    
    axs.plot(gval_ls1, gP_ls1)

    gloc = []
    gloc.append(x_conf)
    gloc.append(x_pval)

    glbl = []
    glbl.append("x: confidence")
    glbl.append("x: p-value")    

    axs.set_xticks(gloc)
    axs.set_xticklabels(glbl)
        
    plt.show()

# peam_show_cdf_chi2(): apresenta o grafico da funcao distribuicao acumulada da Distribuicao Chi-Square
# df - parametro "df" da distribuicao 
# scale - escala para os valores de "x"
# x_conf - posicao da marcacao de "x"
# x_pval - posicao da marcacao do "p-value" 
def peam_show_cdf_chi2(df, scale, x_conf, x_pval):
    # VARS DEF
    #
    min_x = 0.01                     # valor minimo da variavel aleatoria
    max_x = 99.99                    # valor maximo da variavel aleatoria
    dx = 0.01                        # valor do incremento da variavel aleatoria
    x = min_x                        # valor da variavel aleatoria
    x_scale = x * scale              # valor da variavel aleatoria
    y = 0.0                          # valor calculado para a funcao F    
    gval_ls1 = []                    # lista de valores medios da faixa
    gP_ls1 = []                      # lista de valores da funcao distribuicao

    while x <= max_x:
        x_scale = x * scale
        y = peam_func_cdf_chi2(x_scale, df)   
        gval_ls1.append(x_scale)
        gP_ls1.append(y)
        x = x + dx
        
    str_titulo = "cdf-CHI2(x,df=" + str(df) + ") no intervalo [0,1]"

    fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True)
    fig.suptitle(str_titulo)
    
    axs.plot(gval_ls1, gP_ls1)

    gloc = []
    gloc.append(x_conf)
    gloc.append(x_pval)

    glbl = []
    glbl.append("x: confidence")
    glbl.append("x: p-value")    

    axs.set_xticks(gloc)
    axs.set_xticklabels(glbl)
        
    plt.show()



# FIM
