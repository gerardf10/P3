import numpy as np
import matplotlib.pyplot as plt

# === Cargar archivos de datos ===
archivo_pitch = "C:\Users\Usuari\Documents\@ETSETB\PAV\P3\scripts\prueba_pitch.txt"
archivo_pot = "prueba_pot.txt"
archivo_r1norm = "prueba_r1norm.txt"
archivo_rmaxnorm = "prueba_rmaxnorm.txt"

# Leer los datos del pitch (columna 0)
datos_pitch = np.loadtxt(archivo_pitch, usecols=0)

# Leer las otras características
potencia = np.loadtxt(archivo_pot)
r1_normalizado = np.loadtxt(archivo_r1norm)
rmax_normalizado = np.loadtxt(archivo_rmaxnorm)

# Asegurarse de que todos los archivos tengan el mismo número de tramas
num_tramas = min(len(datos_pitch), len(potencia), len(r1_normalizado), len(rmax_normalizado))
tramas = np.arange(num_tramas)

# Recortar los datos para que todos tengan el mismo número de tramas
datos_pitch = datos_pitch[:num_tramas]
potencia = potencia[:num_tramas]
r1_normalizado = r1_normalizado[:num_tramas]
rmax_normalizado = rmax_normalizado[:num_tramas]

# === Graficar ===
fig, axs = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

# 1. Gráfico del Pitch
axs[0].plot(tramas, datos_pitch, label='Pitch (Hz)', color='purple')
axs[0].set_ylabel("Pitch [Hz]")
axs[0].set_title("Contorn de pitch estimat utilizant WaveSurfer")
axs[0].grid(True)
axs[0].legend()

# 2. Gráfico de la Potencia
axs[1].plot(tramas, potencia, label='Potència (r[0]) [dB]', color='black')
axs[1].set_ylabel("r[0] [dB]")
axs[1].set_title("Nivell de potència del senyal")
axs[1].grid(True)
axs[1].legend()

# 3. Gráfico de r1_normalizado
axs[2].plot(tramas, r1_normalizado, label='r[1]/r[0]', color='blue')
axs[2].set_ylabel("r1_normalizat")
axs[2].set_title("Autocorrelació normalitzada (r1_normalizat)")
axs[2].grid(True)
axs[2].legend()

# 4. Gráfico de rmax_normalizado
axs[3].plot(tramas, rmax_normalizado, label='r[lag]/r[0]', color='green')
axs[3].set_ylabel("rmax_normalizat")
axs[3].set_title("Autocorrelació al màxim secundari (rmax_normalitzat)")
axs[3].set_xlabel("Número de trama")
axs[3].grid(True)
axs[3].legend()

# Ajustar el espaciado entre subgráficos
plt.tight_layout()
# Mostrar la figura
plt.show()