#!/bin/bash

#SBATCH --partition=roma
#SBATCH --account=supercdms:default
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8192
#SBATCH --time=0-08:00:00
#SBATCH --gpus 0
#SBATCH --output=slurm-photon-%j.out

# Runs photon_nb_backplane_qpdensity.mac (a 1.4 THz photon in the Nb backplane,
# 2 QPs of 2.895 meV) as a SLURM job.
#
# Submit from the G4Macros directory:
#   sbatch photon.sh [run_dir]
# One photon (one event) with the random seed fixed to 7.
#
# Everything the job writes (the macro copy, the QP step file, the hit and
# primary files, the log) goes into its own run directory, because the QP step
# file is always named QuasiparticleStepInformationFile.txt and is written to the
# working directory. Afterwards, in that directory:
#   root -l -b -q '<path>/AnalysisTools/qp_count_vs_time.C("QuasiparticleStepInformationFile.txt")'

# The SCDMS module is not visible in every batch environment; apptainer is run
# below with the full path of the image, so a missing module is not fatal
module load scdms/V05-02 || echo "WARNING: could not load module scdms/V05-02; continuing"
module list

BaseMacro="${SLURM_SUBMIT_DIR}/photon_nb_backplane_qpdensity.mac"
nEvents=1	# one photon per run
seed=7

# Geometry settings written into the macro copy
Al_bool=false			# /g4cmp/Al false = Nb patterned film
Al_z=100			# Al film thickness [nm]
Nb_z=20				# Nb film thickness [nm] (the backplane stays 200 nm)
pAbsProbPolishedWallSi=0.01	# phonon absorption on the polished Si faces
run_dir="${1:-/sdf/home/s/shilin/mycode/photon_run/run_${SLURM_JOB_ID}}"

# The built example; change this (or set G4EXE before sbatch) to your build
: "${G4EXE:=/sdf/home/s/shilin/mycode/TFR/G4CMP/examples/quasiparticle/quasiparticle-install/bin/g4cmpQuasiparticle}"

# G4CMP environment (only set when not already set in your shell)
: "${G4CMPINSTALL:=/sdf/home/s/shilin/mycode/TFR/G4CMP_install/share/G4CMP}"
: "${G4CMPLIB:=/sdf/home/s/shilin/mycode/TFR/G4CMP_install/lib}"
: "${G4LATTICEDATA:=$G4CMPINSTALL/CrystalMaps}"
export G4CMPINSTALL G4CMPLIB G4LATTICEDATA
export LD_LIBRARY_PATH="$G4CMPLIB:$LD_LIBRARY_PATH"

if [ ! -f "$BaseMacro" ]
then
    echo "Cannot find $BaseMacro: submit from the G4Macros directory"
    exit 1
fi
if [ ! -f "$G4EXE" ]
then
    echo "Cannot find the executable $G4EXE"
    echo "Set G4EXE before sbatch (export G4EXE=/path/to/g4cmpQuasiparticle) or edit this script."
    echo "Find your build with:  find /sdf/home/s/shilin/mycode -name g4cmpQuasiparticle -type f"
    exit 1
fi
if ! command -v apptainer > /dev/null
then
    echo "apptainer is not available in this job"
    exit 1
fi

mkdir -p "$run_dir"
mName="$run_dir/Run_photon.mac"

echo "Base macro:  $BaseMacro"
echo "Run macro:   $mName"
echo "Run dir:     $run_dir"
echo "Executable:  $G4EXE"
echo "Events:      $nEvents"
echo "Seed:        $seed"
echo "Al bool:     $Al_bool"
echo "Al thickness: $Al_z nm"
echo "Nb thickness: $Nb_z nm"
echo "pAbsProbPolishedWallSi: $pAbsProbPolishedWallSi"

cp "$BaseMacro" "$mName"
sed -i "s|^/run/beamOn .*|/run/beamOn ${nEvents}|" "$mName"
sed -i "s|^/random/setSeeds .*|/random/setSeeds ${seed} $((seed + 1))|" "$mName"
sed -i "s|^/g4cmp/Al .*|/g4cmp/Al ${Al_bool}|" "$mName"
sed -i "s|^/g4cmp/filmThicknessAl .*|/g4cmp/filmThicknessAl ${Al_z}|" "$mName"
sed -i "s|^/g4cmp/filmThicknessNb .*|/g4cmp/filmThicknessNb ${Nb_z}|" "$mName"
sed -i "s|^/g4cmp/pAbsProbPolishedWallSi .*|/g4cmp/pAbsProbPolishedWallSi ${pAbsProbPolishedWallSi}|" "$mName"

# The output files (photon_primary.txt, photon_hits.txt and the QP step file)
# are written to the working directory
cd "$run_dir" || exit 1

# The container puts its own G4CMP (/opt/cdms/release/lib) ahead of any
# LD_LIBRARY_PATH exported outside, so set the path inside the container, right
# before the program starts, to use the G4CMP install
apptainer exec /sdf/group/supercdms/software/releases/cdmsfull_V05-02.sif \
  env LD_LIBRARY_PATH="$G4CMPLIB:$LD_LIBRARY_PATH" \
  "$G4EXE" "$mName" > "$run_dir/photon.log" 2>&1

echo "Done. Step file: $run_dir/QuasiparticleStepInformationFile.txt"
