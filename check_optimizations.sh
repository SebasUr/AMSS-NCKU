#!/bin/bash
#
# Script de verificación para las optimizaciones OpenMP en AMSS-NCKU
# Ejecutar antes de compilar para verificar que todo está listo
#

echo "========================================"
echo "AMSS-NCKU OpenMP Optimization Checker"
echo "========================================"
echo ""

ERRORS=0
WARNINGS=0

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check 1: Verify OpenMP directives in diff_new.f90
echo "1. Checking OpenMP directives in diff_new.f90..."
DIFF_OMP_COUNT=$(grep -c '$omp parallel do' AMSS_NCKU_source/diff_new.f90)
if [ "$DIFF_OMP_COUNT" -ge 40 ]; then
    echo -e "   ${GREEN}✓${NC} Found $DIFF_OMP_COUNT OpenMP directives (expected ~45)"
else
    echo -e "   ${RED}✗${NC} Found only $DIFF_OMP_COUNT OpenMP directives (expected ~45)"
    ERRORS=$((ERRORS+1))
fi

# Check 2: Verify OpenMP directives in bssn_rhs.f90
echo "2. Checking OpenMP directives in bssn_rhs.f90..."
BSSN_OMP_COUNT=$(grep -c '$omp parallel do' AMSS_NCKU_source/bssn_rhs.f90)
if [ "$BSSN_OMP_COUNT" -ge 2 ]; then
    echo -e "   ${GREEN}✓${NC} Found $BSSN_OMP_COUNT OpenMP directives"
else
    echo -e "   ${YELLOW}⚠${NC} Found only $BSSN_OMP_COUNT OpenMP directives (expected 2+)"
    WARNINGS=$((WARNINGS+1))
fi

# Check 3: Verify compiler flags in makefile.inc
echo "3. Checking compiler optimization flags..."
if grep -q "march=native" AMSS_NCKU_source/makefile.inc && \
   grep -q "ffast-math" AMSS_NCKU_source/makefile.inc && \
   grep -q "fopenmp" AMSS_NCKU_source/makefile.inc; then
    echo -e "   ${GREEN}✓${NC} Compiler flags are optimized"
else
    echo -e "   ${RED}✗${NC} Compiler flags missing optimizations"
    ERRORS=$((ERRORS+1))
fi

# Check 4: Verify Slurm configuration
echo "4. Checking Slurm job configuration..."
if grep -q "ntasks-per-node=16" run_amss_slurm_longjobs.sh && \
   grep -q "cpus-per-task=4" run_amss_slurm_longjobs.sh; then
    echo -e "   ${GREEN}✓${NC} Slurm hybrid MPI+OpenMP configuration correct"
else
    echo -e "   ${RED}✗${NC} Slurm configuration not optimal"
    ERRORS=$((ERRORS+1))
fi

# Check 5: Verify OMP_NUM_THREADS setting
echo "5. Checking OpenMP environment variables..."
if grep -q "OMP_NUM_THREADS=4" run_amss_slurm_longjobs.sh && \
   grep -q "OMP_PROC_BIND=close" run_amss_slurm_longjobs.sh; then
    echo -e "   ${GREEN}✓${NC} OpenMP environment variables configured"
else
    echo -e "   ${YELLOW}⚠${NC} OpenMP environment variables may need adjustment"
    WARNINGS=$((WARNINGS+1))
fi

# Check 6: Verify MPI process count in Python input
echo "6. Checking MPI process count in AMSS_NCKU_Input.py..."
MPI_PROCS=$(grep "MPI_processes" AMSS_NCKU_Input.py | grep -o '[0-9]\+' | head -1)
if [ "$MPI_PROCS" = "32" ]; then
    echo -e "   ${GREEN}✓${NC} MPI processes = $MPI_PROCS (optimal for 2 nodes × 16 ranks)"
elif [ "$MPI_PROCS" = "128" ]; then
    echo -e "   ${YELLOW}⚠${NC} MPI processes = $MPI_PROCS (old value, will cause oversubscription)"
    echo "      Consider changing to 32 in AMSS_NCKU_Input.py"
    WARNINGS=$((WARNINGS+1))
else
    echo -e "   ${YELLOW}⚠${NC} MPI processes = $MPI_PROCS (unexpected value)"
    WARNINGS=$((WARNINGS+1))
fi

# Check 7: Check if backup files exist
echo "7. Checking for backup files..."
if [ -f "AMSS_NCKU_source/diff_new.f90.bak" ]; then
    echo -e "   ${GREEN}✓${NC} Backup file exists (diff_new.f90.bak)"
else
    echo -e "   ${YELLOW}⚠${NC} No backup file found (consider creating one)"
    WARNINGS=$((WARNINGS+1))
fi

# Check 8: Verify modules are available
echo "8. Checking required modules..."
if command -v module &> /dev/null; then
    if module avail gnu14 2>&1 | grep -q gnu14; then
        echo -e "   ${GREEN}✓${NC} GNU compiler module available"
    else
        echo -e "   ${YELLOW}⚠${NC} GNU compiler module may not be available"
        WARNINGS=$((WARNINGS+1))
    fi
else
    echo -e "   ${YELLOW}⚠${NC} Module system not available (may be OK if not using modules)"
    WARNINGS=$((WARNINGS+1))
fi

echo ""
echo "========================================"
echo "Summary:"
echo "========================================"

if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
    echo -e "${GREEN}✓ All checks passed! Ready to compile and run.${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. cd AMSS_NCKU_source && make clean && make ABE"
    echo "  2. cd .. && sbatch run_amss_slurm_longjobs.sh"
    exit 0
elif [ $ERRORS -eq 0 ]; then
    echo -e "${YELLOW}⚠ $WARNINGS warning(s) found. Review above.${NC}"
    echo "   Code should work but may not be optimal."
    exit 0
else
    echo -e "${RED}✗ $ERRORS error(s) and $WARNINGS warning(s) found.${NC}"
    echo "   Please fix errors before compiling."
    exit 1
fi
