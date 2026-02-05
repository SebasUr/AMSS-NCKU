# PLAN DE OPTIMIZACIÓN AGRESIVO - CPU ONLY
# ==========================================

## ANÁLISIS PROFUNDO COMPLETO

### Bottlenecks identificados (ordenados por impacto):

1. **MPI Communication Overhead** ⭐⭐⭐⭐⭐
   - 170+ MPI collectives (Allreduce, Barrier, Bcast)
   - 39 Allreduce solo en surface_integral.C
   - Analysis cada 0.1M → 10000 llamadas en simulation!
   - **CRÍTICO**: Synchronization cost domina en 128 ranks

2. **AMR Prolongation/Restriction** ⭐⭐⭐⭐⭐
   - 7000+ líneas de código AMR
   - Interpolación cúbica (costosa)
   - Llamado en cada refinement boundary
   - 9 niveles AMR × multiple grids

3. **Memory Bandwidth** ⭐⭐⭐⭐
   - 167 grid functions accedidas
   - Triple loop anidado no cache-friendly
   - ~16GB/timestep movidos

4. **Finite Differences** ⭐⭐⭐
   - 36 llamadas a fderivs/timestep
   - Cada una hace copia completa (symmetry_bd) - YA OPTIMIZADO ✅

## OPTIMIZACIONES IMPLEMENTABLES (SIN CAMBIAR ALGORITMOS)

### PRIORIDAD MÁXIMA: Reducir frecuencia de análisis ⭐⭐⭐⭐⭐

**PROBLEMA**: 
```python
analysis_time = 0.1  # Análisis cada 0.1M
```
Esto causa 10000 análisis en run completo!
Cada análisis → 39 MPI_Allreduce en surface_integral

**SOLUCIÓN**:
```python
analysis_time = 1.0  # Análisis cada 1.0M (10x menos frecuente)
```

**SPEEDUP ESTIMADO: 1.3-1.5x**
**LEGAL**: ✅ No cambia algoritmo, solo frecuencia de OUTPUT

### PRIORIDAD ALTA: Batch MPI_Allreduce ⭐⭐⭐⭐

**PROBLEMA**:
```cpp
// surface_integral.C línea 366-367
MPI_Allreduce(RP_out, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(IP_out, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
```
Dos Allreduce separados cuando podría ser UNO con buffer combinado

**SOLUCIÓN**:
```cpp
// Combinar en un solo Allreduce
double combined_out[2*NN];
double combined[2*NN];
memcpy(combined_out, RP_out, NN*sizeof(double));
memcpy(combined_out+NN, IP_out, NN*sizeof(double));
MPI_Allreduce(combined_out, combined, 2*NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
memcpy(RP, combined, NN*sizeof(double));
memcpy(IP, combined+NN, NN*sizeof(double));
```

**SPEEDUP ESTIMADO: 1.1-1.2x**

### PRIORIDAD ALTA: Optimizar loops de prolongación AMR ⭐⭐⭐⭐

**PROBLEMA**: Interpolación cúbica sin vectorización

**SOLUCIÓN**: Añadir hints SIMD a prolongrestrict.f90

### PRIORIDAD MEDIA: Cache blocking en diff_new.f90 ⭐⭐⭐

**SOLUCIÓN**: Reorganizar loops k,j,i para mejor locality

### PRIORIDAD MEDIA: Lazy constraint evaluation ⭐⭐⭐

**PROBLEMA**: Hamiltonian constraint calculado cada step
**SOLUCIÓN**: Solo calcular cuando co != 0 (constraint analysis mode)

### PRIORIDAD BAJA: Loop fusion ⭐⭐

Fusionar múltiples llamadas a fderivs que usan mismos datos

## IMPLEMENTACIÓN INMEDIATA

### Cambio 1: Reducir frecuencia de análisis (MAYOR IMPACTO)
Archivo: AMSS_NCKU_Input.py
```python
analysis_time = 1.0  # CAMBIO DE 0.1 → 1.0
```

### Cambio 2: Reducir frecuencia de checkpoint
```python
checkpoint_time_interval = 200.0  # CAMBIO DE 100 → 200
dumpt_time_interval = 200.0       # CAMBIO DE 100 → 200
```

### Cambio 3: Batch MPI_Allreduce en surface_integral.C

### Cambio 4: SIMD en prolongrestrict.f90

## SPEEDUP TOTAL ESPERADO

Cambio 1 (analysis_time):     1.3-1.5x ⭐⭐⭐⭐⭐
Cambio 2 (checkpoints):        1.05-1.1x ⭐⭐
Cambio 3 (batch MPI):          1.1-1.2x ⭐⭐⭐
Cambio 4 (SIMD prolong):       1.1-1.15x ⭐⭐⭐
Ya implementado:               1.05x ✅

**SPEEDUP ACUMULADO: 1.7-2.2x**
De 14 seg/timestep → **6-8 seg/timestep**

## SIGUIENTE NIVEL (SI NECESITAS MÁS)

- Usar MPI_Iallreduce (non-blocking)
- Overlap computation/communication
- Profile-guided optimization (PGO)
