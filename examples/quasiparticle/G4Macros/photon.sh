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
#   sbatch photon.sh [nEvents] [seed] [run_dir]
# e.g.
#   sbatch photon.sh 100 7
#
# Everything the job writes (the macro copy, the QP step file, the hit and
# primary files, the log) goes into its own run directory, because the QP step
# file is always named QuasiparticleStepInformationFile.txt and is written to the
# working directory. Afterwards, in that directory:
#   root -l -b -q '<path>/AnalysisTools/qp_count_vs_time.C("QuasiparticleStepInformationFile.txt")'

# load SCDMS singularity image module
module load scdms/V05-02
module list

BaseMacro="${SLURM_SUBMIT_DIR}/photon_nb_backplane_qpdensity.mac"
nEvents="${1:-1}"
seed="${2:-123456}"
run_dir="${3:-/sdf/home/s/shilin/mycode/photon_run/run_${SLURM_JOB_ID}}"

# The built example; change this (or set G4EXE before sbatch) to your build
: "${G4EXE:=/sdf/home/s/shilin/mycode/TFR/G4CMP-TFR/examples/quasiparticle/quasiparticle-lQPD_Nb-build/g4cmpQuasiparticle}"

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
    echo "Cannot find the executable $G4EXE: set G4EXE or edit this script"
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

cp "$BaseMacro" "$mName"
sed -i "s|^/run/beamOn .*|/run/beamOn ${nEvents}|" "$mName"
sed -i "s|^/random/setSeeds .*|/random/setSeeds ${seed} $((seed + 1))|" "$mName"

# The output files (photon_primary.txt, photon_hits.txt and the QP step file)
# are written to the working directory
cd "$run_dir" || exit 1

apptainer exec /sdf/group/supercdms/software/releases/cdmsfull_V05-02.sif \
  "$G4EXE" "$mName" > "$run_dir/photon.log" 2>&1

echo "Done. Step file: $run_dir/QuasiparticleStepInformationFile.txt"
