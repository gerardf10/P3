import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wav

# === Configuració ===
arxiu_audio = 'prueba.wav'
duracio_quadro = 0.030  # 30 ms
temps_inici_quadro = 0.770  # Moment on comença una part de la senyal (obtingut amb Wavesurfer)

# === Carregar el fitxer WAV ===
taxa_muestreo, dades = wav.read(arxiu_audio)
if dades.ndim > 1:
    dades = dades[:, 0]  # Convertim a mono si és estèreo

# === Extraure el segment de 30 ms ===
muestra_inici = int(temps_inici_quadro * taxa_muestreo)
longitud_quadro = int(duracio_quadro * taxa_muestreo)
quadro = dades[muestra_inici:muestra_inici + longitud_quadro]
quadro = quadro.astype(np.float32)
quadro = quadro - np.mean(quadro)  # Eliminem el component DC

# === Implementació de l'autocorrelació ===
def autocorrelacio(x):
    N = len(x)
    r = np.zeros(N)
    
    # Per a cada retard l (lag)
    for l in range(N):
        suma = 0.0
        # r_{xx}[m] = (1/N) * sum_0^{N-l} x[n] * x[n+l]
        # Recórrer la senyal fins a l'índex N - l per evitar accedir fora de rang
        for n in range(N - l):
            suma += x[n] * x[n + l]
        
        # Autocorrelació: es divideix pel nombre total de mostres N
        r[l] = suma / N

    # Ajustar r[0] per evitar problemes posteriors (com divisió per zero o logaritmes)
    if r[0] == 0.0:
        r[0] = 1e-10
    return r

# === Calcular l'autocorrelació ===
autocorrelacio = autocorrelacio(quadro)

# === Identificar el màxim secundari (lag != 0) ===
autocorrelacio[0] = 1e-10  # Evitar detectar el màxim en el temps zero
lag_minim = int(taxa_muestreo / 500)  # Pitch màxim: 500 Hz
lag_maxim = int(taxa_muestreo / 50)   # Pitch mínim: 50 Hz

rang_lags = autocorrelacio[lag_minim:lag_maxim]
lag = np.argmax(rang_lags) + lag_minim
periodo_pitch = lag
frequencia_pitch = taxa_muestreo / periodo_pitch

# === Graficar ===
temps = np.arange(longitud_quadro) / taxa_muestreo * 1000  # Temps en ms
lags = np.arange(len(autocorrelacio))

plt.figure(figsize=(10, 6))

# Subgraf 1: senyal temporal
plt.subplot(2, 1, 1)
plt.plot(temps, quadro)
plt.title(f"Senyal temporal (30 ms) extreta de {arxiu_audio} — Pitch = {frequencia_pitch:.2f} Hz")
plt.xlabel("Temps [ms]")
plt.ylabel("Amplitud")
plt.axvline(x=periodo_pitch / taxa_muestreo * 1000, color='r', linestyle='--', label='Període del pitch')
plt.legend()

# Subgraf 2: autocorrelació
plt.subplot(2, 1, 2)
plt.plot(lags, autocorrelacio)
plt.title("Autocorrelació")
plt.xlabel("Lag [mostres]")
plt.ylabel("r[lag]")
plt.axvline(x=lag, color='r', linestyle='--', label='Primer màxim secundari')
plt.legend()

plt.tight_layout()
plt.show()