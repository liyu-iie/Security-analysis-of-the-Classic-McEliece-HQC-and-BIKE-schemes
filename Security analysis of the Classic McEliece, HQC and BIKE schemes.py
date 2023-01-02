import math

def H(x):
    h = -x * math.log2(x) - (1-x) * math.log2(1-x)
    return h

def COM(n, w):
    x = w/n
    h = -x * math.log2(x) - (1-x) * math.log2(1-x)
    y = n*h
    return y

# McEliece_level1 = {"name": "McEliece C1", "n": 3488, "k": 2720, "w": 64}
# McEliece_level3 = {"name": "McEliece C3", "n": 4608, "k": 3360, "w": 96}
# McEliece_level5a = {"name": "McEliece C5a", "n": 6688, "k": 5024, "w": 128}
# McEliece_level5b = {"name": "McEliece C5b", "n": 6960, "k": 5413, "w": 119}
# McEliece_level5c = {"name": "McEliece C5c", "n": 8192, "k": 6528, "w": 128}

# BIKE_level1 = {"name": "BIKE C1", "n": 24646, "k": 12323, "w": 134, "w_k": 142}
# BIKE_level3 = {"name": "BIKE C3", "n": 49318, "k": 24659, "w": 199, "w_k": 206}
# BIKE_level5 = {"name": "BIKE C5", "n": 81946, "k": 40973, "w": 264, "w_k": 274}

# HQC_level1 = {"name": "HQC C1", "n": 35338, "k": 17669, "w": 132, "w_e": 75}
# HQC_level3 = {"name": "HQC C3", "n": 71702, "k": 35851, "w": 200, "w_e": 114}
# HQC_level5 = {"name": "HQC C5", "n": 115274, "k": 57637, "w": 262, "w_e": 149}

####==============================
#### ========== 4-list ===========

####==============================
### McEliece 128-bit
n = 3488
k = 2720
t = 64
DOOM = 1
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
w = 3
v = 5
l = 57
## --- Cube-root memory penalty
# w = 2
# v = 13
# l = 49

####==============================
### McEliece 192-bit
# n = 4608
# k = 3360
# t = 96
# DOOM = 1
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 3
# v = 6
# l = 60
## --- Cube-root memory penalty
# w = 2
# v = 14
# l = 52

####==============================
### McEliece 256-bit
# n = 6688
# k = 5024
# t = 128
# DOOM = 1
## --- M<60, logarithmic memory penalty and cube-root memory penalty
# w = 3
# v = 5
# l = 63

####==============================
### McEliece 256-bit
# n = 6960
# k = 5413
# t = 119
# DOOM = 1
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# v = 0
# l = 76
## --- Cube-root memory penalty
# w = 3
# v = 5
# l = 63

####==============================
### McEliece 256-bit
# n = 8192
# k = 6528
# t = 128
# DOOM = 1
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# v = 0
# l = 78
## --- Cube-root memory penalty
# w = 3
# v = 5
# l = 65

####==============================
### HQC 128-bit --- M<60, logarithmic memory penalty and cube-root memory penalty
# n = 35338
# k = 17669
# t = 132
# DOOM = math.sqrt(k)
# w = 2
# v = 20
# l = 68

####==============================
### HQC 192-bit --- M<60, logarithmic memory penalty and cube-root memory penalty
# n = 71702
# k = 35851
# t = 200
# DOOM = math.sqrt(k)
# w = 2
# v = 21
# l = 73

####==============================
### HQC 256-bit --- M<60, logarithmic memory penalty and cube-root memory penalty
# n = 115274
# k = 57637
# t = 262
# DOOM = math.sqrt(k)
# w = 2
# v = 22
# l = 76

####==============================
### BIKE 128-bit, "w_k": 142
# n = 24646
# k = 12323
# w = 2
# v = 19
# l = 65
# DOOM = math.sqrt(k)
## --- message security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 134
## --- ley security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 142

####==============================
### BIKE 192-bit, "w_k": 206
# n = 49318
# k = 24659
# w = 2
# v = 20
# l = 70
# DOOM = math.sqrt(k)
## --- message security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 199
## --- ley security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 206

####==============================
### BIKe 256-bit, "w_k": 274
# n = 81946
# k = 40973
# w = 2
# v = 21
# l = 73
# DOOM = math.sqrt(k)
## --- message security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 264
## --- ley security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 274

###----------------Theorem 2 in paper --------------------------------------------------
K = 4
m = math.ceil(math.log2(math.comb(int((k + l) / K), w)))

Gauss = (n - k) * (n - k) * (n + 1)

P11 = min(math.comb(n, t) // pow(2, l), pow(2, n - k - l ))
P12 = max(P11, 1)
P1 = math.comb(n - k - l, t - K * w) / P12

# the exhaustive A4 algorithm: l=2m+v
Num = max(1 / (P1 * pow(2, m)), 1)

Iter = max(Gauss, pow(2, v + m))

T = Num * Iter / DOOM  ### DOOM=1 for McEliece; DOOM=math.sqrt(k) for HQC and BIKE

logT = round(math.log2(T), 2)

M = max(math.ceil( math.log2(n * (n - k))), m)


print("l=", l, "  m=", m, "  p=", 4*w, "  v=", v, "  M=", M,"  T=", logT, "  T1=", logT + math.log2(M), "  T2=", logT + M/3)



####========================================================================================

#### ========== 2-list ===========

####==============================
### McEliece 128-bit
# n = 3488
# k = 2720
# t = 64
# DOOM = 1
# w = 4
# l = 38

####==============================
### McEliece 192-bit
# n = 4608
# k = 3360
# t = 96
# DOOM = 1
# w = 4
# l = 39

####==============================
### McEliece 256-bit
# n = 6688
# k = 5024
# t = 128
# DOOM = 1
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 6
# l = 59
## --- Cube-root memory penalty
# w = 3
# l = 32

####==============================
### McEliece 256-bit
# n = 6960
# k = 5413
# t = 119
# DOOM = 1
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 5
# l = 51
## --- Cube-root memory penalty
# w = 3
# l = 32

####==============================
### McEliece 256-bit
# n = 8192
# k = 6528
# t = 128
# DOOM = 1
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 5
# l = 52
## --- Cube-root memory penalty
# w = 6
# l = 61
## --- Cube-root memory penalty
# w = 3
# l = 33

####==============================
### HQC 128-bit --- M<60, logarithmic memory penalty and cube-root memory penalty
# n = 35338
# k = 17669
# t = 132
# DOOM = math.sqrt(k)
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 48
## --- Cube-root memory penalty
# w = 3
# l = 37

####==============================
### HQC 192-bit --- M<60, logarithmic memory penalty and cube-root memory penalty
# n = 71702
# k = 35851
# t = 200
# DOOM = math.sqrt(k)
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 52
## --- Cube-root memory penalty
# w = 3
# l = 40

####==============================
### HQC 256-bit --- M<60, logarithmic memory penalty and cube-root memory penalty
# n = 115274
# k = 57637
# t = 262
# DOOM = math.sqrt(k)
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 55
## --- Cube-root memory penalty
# w = 3
# l = 42

####==============================
### BIKE 128-bit, "w_k": 142
# n = 24646
# k = 12323
# DOOM = math.sqrt(k)

## --- message security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 134
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 46
## --- Cube-root memory penalty
# w = 3
# l = 36

## --- ley security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 142
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 46
## --- Cube-root memory penalty
# w = 3
# l = 36

####==============================
### BIKE 192-bit, "w_k": 206
# n = 49318
# k = 24659
# DOOM = math.sqrt(k)

## --- message security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 199
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 50
## --- Cube-root memory penalty
# w = 3
# l = 39

## --- ley security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 206
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 50
## --- Cube-root memory penalty
# w = 3
# l = 39

####==============================
### BIKe 256-bit, "w_k": 274
# n = 81946
# k = 40973
# DOOM = math.sqrt(k)

## --- message security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 264
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 53
## --- Cube-root memory penalty
# w = 3
# l = 41

## --- ley security --- M<60, logarithmic memory penalty and cube-root memory penalty
# t = 274
## --- No memory penalty i.e., memory M<60 and logarithmic memory penalty
# w = 4
# l = 53
## --- Cube-root memory penalty
# w = 3
# l = 41

###----------------Theorem 2 in paper --------------------------------------------------
# K = 2
# m = math.ceil(math.log2(math.comb(int((k + l) / K), w)))
#
# Gauss = (n - k) * (n - k) * (n + 1)
#
# P11 = min(math.comb(n, t) // pow(2, l), pow(2, n - k - l ))
# P12 = max(P11, 1)
# P1 = math.comb(n - k - l, t - K * w) / P12
#
# # the exhaustive A4 algorithm: l=2m+v
# Num = max(1 / (P1 * pow(2, m)), 1)
#
# Iter = max(Gauss, pow(2, m))
#
# T = Num * Iter / DOOM  ### DOOM=1 for McEliece; DOOM=math.sqrt(k) for HQC and BIKE
#
# logT = round(math.log2(T), 2)
#
# M = max(math.ceil( math.log2(n * (n - k))), m)
#
#
# print("l=", l, "  m=", m, "  p=", 2*w, "  M=", M,"  T=", logT, "  T1=", logT + math.log2(M), "  T2=", logT + M/3)
