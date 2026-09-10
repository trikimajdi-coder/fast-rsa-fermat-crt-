import ctypes
import os
import math
import time

lib_path = os.path.abspath("./libcrt_filter_31_omp.so")
if not os.path.exists(lib_path):
    raise FileNotFoundError(f"Fichier {lib_path} introuvable. Compilez-le avec 'make'.")

lib = ctypes.CDLL(lib_path)

class CRTResult(ctypes.Structure):
    _fields_ = [
        ("iterations", ctypes.c_ulonglong),
        ("execution_time", ctypes.c_double),
        ("success", ctypes.c_int)
    ]

lib.fast_crt_factorize_31_omp.argtypes = [
    ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_int), ctypes.c_int,
    ctypes.POINTER(ctypes.c_ulonglong),
    ctypes.c_ulonglong,
    ctypes.c_int,
    ctypes.c_char_p
]
lib.fast_crt_factorize_31_omp.restype = CRTResult

def solve_rsa_crt_omp(N_str: str, k_max: int = 5000):
    N = int(N_str)
    modules = [6, 35, 143, 17, 19, 23, 29, 31]
    M = math.prod(modules)
    
    valid_res = []
    for m in modules:
        quad_res = { (y * y) % m for y in range(m) }
        res_m = [rx for rx in range(m) if ((rx * rx - N) % m) in quad_res]
        valid_res.append(res_m)
        
    W = []
    for m_i in modules:
        C_i = M // m_i
        W.append((C_i * pow(C_i, -1, m_i)) % M)

    c_r = [(ctypes.c_int * len(vr))(*vr) for vr in valid_res]
    c_W = (ctypes.c_ulonglong * len(W))(*W)
    out_x_buffer = ctypes.create_string_buffer(2048)

    res = lib.fast_crt_factorize_31_omp(
        N_str.encode('utf-8'),
        c_r[0], len(valid_res[0]),
        c_r[1], len(valid_res[1]),
        c_r[2], len(valid_res[2]),
        c_r[3], len(valid_res[3]),
        c_r[4], len(valid_res[4]),
        c_r[5], len(valid_res[5]),
        c_r[6], len(valid_res[6]),
        c_r[7], len(valid_res[7]),
        c_W,
        M,
        k_max,
        out_x_buffer
    )

    if res.success:
        x_found = int(out_x_buffer.value.decode('utf-8'))
        y_found = math.isqrt(x_found**2 - N)
        return {
            "success": True,
            "x": x_found,
            "p": x_found - y_found,
            "q": x_found + y_found,
            "execution_time_sec": res.execution_time,
            "iterations": res.iterations
        }
    return {"success": False, "execution_time_sec": res.execution_time, "iterations": res.iterations}

if __name__ == "__main__":
    p_test = 1000000007
    q_test = 1000000009
    N_test = str(p_test * q_test)
    
    t0 = time.perf_counter()
    resultat = solve_rsa_crt_omp(N_test, k_max=10000)
    t1 = time.perf_counter()

    if resultat["success"]:
        print(f"[SUCCÈS] P={resultat['p']}, Q={resultat['q']} en {resultat['execution_time_sec']*1000:.3f} ms")
    else:
        print("[ÉCHEC]")
