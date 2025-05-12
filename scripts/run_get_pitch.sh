#!/bin/bash

# Establecemos que el código de retorno de un pipeline sea el del último programa con código de retorno
# distinto de cero, o cero si todos devuelven cero.
set -o pipefail

# override to the Linux‐side binary
GETF0="/home/guillem/PAV/bin/get_pitch"

# Windows‐side data directory (this still works for reading .wav/.f0ref)
WAV_DIR="/home/guillem/PAV/p3/pitch_db/train"

pitch_evaluate="/home/guillem/PAV/bin/pitch_evaluate"

# Per a cada fitxer .wav en el directori
for fwav in "$WAV_DIR"/*.wav; do
    ff0=${fwav/.wav/.f0}
    echo "$GETF0 $fwav $ff0 ----"
    $GETF0 --clip --median --lowpass=350 $fwav $ff0 > /dev/null || { echo -e "\nError in $GETF0 $fwav $ff0" && exit 1; }
done

$pitch_evaluate "$WAV_DIR"/*.f0ref

exit 0