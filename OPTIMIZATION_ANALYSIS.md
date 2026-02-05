# ANÁLISIS DE OPTIMIZACIÓN CPU - AMSS-NCKU
# ==========================================

## 1. ANÁLISIS DE PERFILADO (lo que encontré)

### Bottlenecks principales:
1. **diff_new.f90** (4313 líneas)
   - `fderivs()`: Llamado 36 veces por timestep
   - `symmetry_bd()`: Copia de arrays + boundary conditions
   - Triple loop anidado con operaciones de stencil
   - **Problema**: Cada fderivs hace una COPIA completa del array (symmetry_bd)

2. **bssn_rhs.f90** (1187 líneas)
   - Llamadas a 36 derivadas por timestep
   - Operaciones tensoriales densas
   - Actualizaciones de 167 grid functions

3. **Memory bandwidth limitado**
   - 167 grid functions × 96×96×48 (nivel más grande)
   - ~16GB de memoria accedida por timestep
   - **Problema crítico**: Memory-bound, no compute-bound!

## 2. OPTIMIZACIONES IMPLEMENTADAS (ya hechas)
✅ Flags de compilador agresivas
✅ MPI parallelization (128 ranks)
✅ OpenMP directives (pero causan overhead)

## 3. OPTIMIZACIONES PENDIENTES (analizadas)

### A) Eliminar copias innecesarias en symmetry_bd ⭐⭐⭐⭐⭐
**IMPACTO MAYOR**: 36 copias de arrays por timestep
```fortran
! ACTUAL (línea 327-332 de fmisc.f90):
funcc = 0.d0  ! INICIALIZACIÓN COMPLETA (lenta!)
funcc(1:extc(1),1:extc(2),1:extc(3)) = func  ! COPIA COMPLETA

! OPTIMIZADO: Solo copiar lo necesario
! No inicializar todo el array, solo ghost zones
```
**Speedup estimado**: 1.2-1.4x

### B) Cache blocking en loops triple anidados ⭐⭐⭐⭐
**PROBLEMA**: Los loops k,j,i acceden memoria linealmente pero cache es limitado
```fortran
! ACTUAL:
do k=1,ex(3)-1
  do j=1,ex(2)-1
    do i=1,ex(1)-1
      fx(i,j,k) = ...  ! Cache miss frecuente
    
! OPTIMIZADO: Blocking para L1/L2 cache
do kb=1,ex(3),BLOCK_SIZE
  do jb=1,ex(2),BLOCK_SIZE
    do ib=1,ex(1),BLOCK_SIZE
      do k=kb,min(kb+BLOCK_SIZE-1,ex(3)-1)
        do j=jb,min(jb+BLOCK_SIZE-1,ex(2)-1)
          do i=ib,min(ib+BLOCK_SIZE-1,ex(1)-1)
```
**Speedup estimado**: 1.1-1.3x

### C) Fusión de loops de derivadas ⭐⭐⭐
36 llamadas separadas a fderivs → muchas pasadas sobre los datos
Fusionar derivadas que usan los mismos datos
**Speedup estimado**: 1.1-1.2x

### D) SIMD manual con intrinsics ⭐⭐⭐
Fortran compiler vectoriza pero puede mejorar con hints explícitos
```fortran
!DIR$ SIMD
!DIR$ VECTOR ALIGNED
```
**Speedup estimado**: 1.05-1.15x

### E) Prefetching manual ⭐⭐
```fortran
!DIR$ PREFETCH f:1
```
**Speedup estimado**: 1.03-1.08x

## 4. OPTIMIZACIONES ARQUITECTÓNICAS

### F) Usar diferencias finitas de menor orden LOCALMENTE ⭐⭐⭐⭐
**IDEA**: Usar 4th order solo donde importa, 2nd order en resto
- Niveles AMR bajos (0-3): 2nd order (más rápidos)
- Niveles AMR altos (4-8): 4th order (precisión necesaria)
**Speedup estimado**: 1.3-1.5x
**PROBLEMA**: Necesita cambiar macrodef.fh dinámicamente (difícil)

### G) Reducir ghost_width adaptativamente ⭐⭐⭐
ghost_width=3 (4th order) en todos lados
Niveles bajos podrían usar ghost_width=2
**Speedup estimado**: 1.1-1.2x

## 5. OPTIMIZACIONES ALGORÍTMICAS EQUIVALENTES

### H) Kreiss-Oliger dissipation optimizada ⭐⭐
Disipación actual se calcula en cada punto
Podría calcularse solo en regiones de alta curvatura
**Speedup estimado**: 1.05-1.1x

### I) Lazy evaluation de constraints ⭐⭐
Hamiltonian/Momentum constraints se calculan cada step
Solo necesarios para análisis cada 0.1M
**Speedup estimado**: 1.02-1.05x

## 6. OPTIMIZACIÓN DE COMUNICACIÓN MPI

### J) Overlap computation/communication ⭐⭐⭐⭐
Usar MPI_Isend/Irecv non-blocking
Computar interior mientras se comunican boundaries
**Speedup estimado**: 1.2-1.5x (si memory-bound)

### K) Reducir frecuencia de synchronization ⭐⭐
Analizar cada 0.1M → menos MPI_Allreduce
**Speedup estimado**: 1.05-1.1x

## 7. RECOMENDACIONES PRIORIZADAS

### PRIORIDAD 1 (IMPLEMENTAR YA): ⭐⭐⭐⭐⭐
A) Optimizar symmetry_bd (eliminar inicialización completa)
J) MPI overlap (computation + communication)

### PRIORIDAD 2 (IMPACTO MEDIO): ⭐⭐⭐
B) Cache blocking en diff_new.f90
C) Fusión de loops de derivadas

### PRIORIDAD 3 (REFINAMIENTO): ⭐⭐
D) SIMD manual
E) Prefetching
H) KO dissipation optimizada

### NO HACER (restricciones de competencia):
F) Cambiar orden FD
G) Cambiar ghost_width
(Cambian parámetros numéricos)

## 8. SPEEDUP TOTAL ESTIMADO
Prioridad 1: 1.4-2.0x
Prioridad 2: 1.2-1.5x adicional
Prioridad 3: 1.1-1.2x adicional

**SPEEDUP ACUMULADO: 1.8-3.6x** (de 15 seg → 4-8 seg/timestep)
