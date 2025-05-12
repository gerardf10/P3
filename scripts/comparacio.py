import numpy as np
import matplotlib.pyplot as plt

# === Cargar las estimaciones de pitch de ambos métodos ===
archivo_pitch_propio = "prueba.f0.txt"
archivo_pitch_ws = "prueba_pitch.txt"

# Leer los datos de ambos archivos
pitch_propio = np.loadtxt(archivo_pitch_propio)
pitch_ws = np.loadtxt(archivo_pitch_ws, usecols=0)

# Asegurarse de que ambos tengan el mismo número de tramas
num_tramas = min(len(pitch_propio), len(pitch_ws))
indices_tramas = np.arange(num_tramas)

# Recortar los datos para que tengan el mismo tamaño
pitch_propio = pitch_propio[:num_tramas]
pitch_ws = pitch_ws[:num_tramas]

# === Crear los gráficos ===
plt.figure(figsize=(12, 6))
plt.plot(indices_tramas, pitch_propio, label="Estimació pròpia (get_pitch)", color='blue')
plt.plot(indices_tramas, pitch_ws, label="Estimació WaveSurfer", color='orange', linestyle='--')
plt.title("Comparativa entre l'estimador 'get_pitch' i WaveSurfer")
plt.xlabel("Índex de trama")
plt.ylabel("Pitch [Hz]")
plt.legend(loc='upper right')
plt.grid(True)
plt.tight_layout()
plt.show()
