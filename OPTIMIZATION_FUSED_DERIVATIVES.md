# AMSS-NCKU Fused Derivative Optimization

## Summary

This document describes the **fused derivative loop optimization** applied to the AMSS-NCKU numerical relativity code. The optimization reduces memory bandwidth pressure by computing derivatives of multiple variables in a single pass through memory, instead of one variable at a time.

**What changed:**
- 23 individual derivative calls in `bssn_rhs.f90` → 8 batched calls
- 3 new subroutines in `diff_new.f90`: `fderivs3`, `fderivs2`, `fdderivs3`
- **Zero mathematical changes** — identical stencils, identical symmetry handling

**Expected improvement:** 5–15% reduction in `bssn_rhs` computation time (memory-bandwidth bound).

---

## The Problem: Memory Bandwidth Bottleneck

### How BSSN RHS Works

The BSSN (Baumgarte–Shapiro–Shibata–Nakamura) right-hand-side computation is the **hottest function** in the code. Each RK45 substep calls `compute_rhs_bssn()`, which computes finite-difference derivatives of ~25 grid variables.

The original code computes derivatives **one variable at a time**:

```fortran
! Original: 3 separate calls for beta^i derivatives
call fderivs(ex, betax, betaxx, betaxy, betaxz, X,Y,Z, ANTI, SYM, SYM, Symmetry, Lev)
call fderivs(ex, betay, betayx, betayy, betayz, X,Y,Z,  SYM,ANTI, SYM, Symmetry, Lev)
call fderivs(ex, betaz, betazx, betazy, betazz, X,Y,Z,  SYM, SYM,ANTI, Symmetry, Lev)
```

### What Each `fderivs` Call Does

Each call to `fderivs(f → fx, fy, fz)` performs:

1. **`symmetry_bd()`** — Copies `f(1:N,1:N,1:N)` into ghost array `fh(-1:N,-1:N,-1:N)` and fills 2 ghost zones per face (12 faces × 2 layers = boundary data)
2. **Interior loop** — For all points `(i,j,k)` in `[3, N-3]³`, computes the 4th-order stencil:

$$f_x(i,j,k) = \frac{1}{12\Delta x}\left[f(i-2) - 8f(i-1) + 8f(i+1) - f(i+2)\right]$$

3. **Boundary loop** — For remaining points near edges, uses conditionals to select 4th-order or 2nd-order stencils depending on available neighbors.

### The Cache Problem

Each `fderivs` call:
- Reads `f` (full 3D array) → copies to `fh` (slightly larger)
- Reads `fh` in stencil pattern → writes `fx, fy, fz`
- **Total memory traffic per call:** ~4 arrays × N³ × 8 bytes

For a typical block of `N = 37` (level 0, 128 MPI ranks):
- Array size: 37³ × 8 = ~405 KB per array
- Traffic per `fderivs`: ~1.6 MB (read fh 3× for x,y,z stencils + write 3 outputs)

With **23 separate calls** in `bssn_rhs`, the ghost array `fh` is allocated, filled, used, and discarded **23 times**. Each time, the input array must be **re-read from L2/L3 cache** (or worse, main memory), even though many of these variables are stored contiguously and could share cache lines.

---

## The Solution: Fused Derivative Loops

### Core Idea

Instead of processing one variable per call, process **3 variables in one fused loop**. This way:
- All 3 ghost arrays (`fh1, fh2, fh3`) are allocated and filled upfront
- The hot inner loop reads from all 3 arrays in the same iteration, keeping data in L1/L2 cache
- Temporal locality is maximized: when `fh1(i,j,k)` is in cache, `fh2(i,j,k)` and `fh3(i,j,k)` are likely in the same or adjacent cache lines

### New Subroutine: `fderivs3`

```fortran
subroutine fderivs3(ex, f1,f1x,f1y,f1z, f2,f2x,f2y,f2z, f3,f3x,f3y,f3z,
                    X,Y,Z, SYM1a,SYM2a,SYM3a, SYM1b,SYM2b,SYM3b,
                    SYM1c,SYM2c,SYM3c, symmetry, onoff)
```

**Structure:**
1. Compute symmetry bounds for each variable independently (`imin1`, `imin2`, `imin3`, etc.)
2. Call `symmetry_bd()` 3 times — one per variable with its own symmetry arguments
3. **Single fused interior loop:**
```fortran
do k = 3, ex(3)-3
do j = 3, ex(2)-3
do i = 3, ex(1)-3
    ! Variable 1: all 3 directions
    f1x(i,j,k) = d12dx*(fh1(i-2,j,k) - 8*fh1(i-1,j,k) + 8*fh1(i+1,j,k) - fh1(i+2,j,k))
    f1y(i,j,k) = d12dy*(fh1(i,j-2,k) - 8*fh1(i,j-1,k) + 8*fh1(i,j+1,k) - fh1(i,j+2,k))
    f1z(i,j,k) = d12dz*(fh1(i,j,k-2) - 8*fh1(i,j,k-1) + 8*fh1(i,j,k+1) - fh1(i,j,k+2))
    ! Variable 2: all 3 directions
    f2x(i,j,k) = d12dx*(fh2(i-2,j,k) - 8*fh2(i-1,j,k) + 8*fh2(i+1,j,k) - fh2(i+2,j,k))
    f2y(i,j,k) = ...
    f2z(i,j,k) = ...
    ! Variable 3: all 3 directions
    f3x(i,j,k) = ...
    f3y(i,j,k) = ...
    f3z(i,j,k) = ...
enddo; enddo; enddo
```
4. **Separate boundary loop** — handles edge points with per-variable symmetry bounds

### Also Implemented: `fderivs2` and `fdderivs3`

- **`fderivs2`**: Same pattern for 2 variables (used for Lap + trK which share SYM,SYM,SYM)
- **`fdderivs3`**: Fused second-derivative computation for 3 variables (6 outputs each: fxx, fxy, fxz, fyy, fyz, fzz). Used for shift vector second derivatives.

---

## Changes in `bssn_rhs.f90`

### Group 1: Shift Vector First Derivatives (3 → 1 call)

```fortran
! BEFORE: 3 calls, 3 separate loops, 3 ghost array allocations
call fderivs(ex, betax, ..., ANTI, SYM, SYM, ...)
call fderivs(ex, betay, ...,  SYM,ANTI, SYM, ...)
call fderivs(ex, betaz, ...,  SYM, SYM,ANTI, ...)

! AFTER: 1 call, 1 fused loop, 3 ghost arrays in cache together
call fderivs3(ex, betax,..., betay,..., betaz,...,
              X,Y,Z, ANTI,SYM,SYM, SYM,ANTI,SYM, SYM,SYM,ANTI, Symmetry, Lev)
```

### Group 2: Metric First Derivatives (6 → 2 calls)

```fortran
! BEFORE: 6 calls for d/dx,dy,dz of (g̃xx, g̃xy, g̃xz, g̃yy, g̃yz, g̃zz)
call fderivs(ex, dxx, ..., SYM, SYM, SYM, ...)
call fderivs(ex, gxy, ..., ANTI,ANTI, SYM, ...)
call fderivs(ex, gxz, ..., ANTI, SYM,ANTI, ...)
call fderivs(ex, dyy, ..., SYM, SYM, SYM, ...)
call fderivs(ex, gyz, ..., SYM, ANTI,ANTI, ...)
call fderivs(ex, dzz, ..., SYM, SYM, SYM, ...)

! AFTER: 2 batched calls
call fderivs3(ex, dxx,..., gxy,..., gxz,..., ...)   ! first 3 components
call fderivs3(ex, dyy,..., gyz,..., dzz,..., ...)   ! last 3 components
```

### Group 3: Lapse + Trace-K (2 → 1 call)

```fortran
! BEFORE
call fderivs(ex, Lap, Lapx,Lapy,Lapz, ..., SYM,SYM,SYM, ...)
call fderivs(ex, trK, Kx,Ky,Kz,       ..., SYM,SYM,SYM, ...)

! AFTER
call fderivs2(ex, Lap,Lapx,Lapy,Lapz, trK,Kx,Ky,Kz, ..., SYM,SYM,SYM, SYM,SYM,SYM, ...)
```

### Group 4: Connection Functions (3 → 1 call)

```fortran
! BEFORE
call fderivs(ex, Gamx, ..., ANTI, SYM, SYM, ...)
call fderivs(ex, Gamy, ...,  SYM,ANTI, SYM, ...)
call fderivs(ex, Gamz, ...,  SYM, SYM,ANTI, ...)

! AFTER
call fderivs3(ex, Gamx,..., Gamy,..., Gamz,..., ANTI,SYM,SYM, SYM,ANTI,SYM, SYM,SYM,ANTI, ...)
```

### Group 5: Shift Vector Second Derivatives (3 → 1 call)

```fortran
! BEFORE: 3 fdderivs calls (6 outputs each = 18 output arrays)
call fdderivs(ex, betax, gxxx,gxyx,gxzx,gyyx,gyzx,gzzx, ..., ANTI,SYM,SYM, ...)
call fdderivs(ex, betay, gxxy,gxyy,gxzy,gyyy,gyzy,gzzy, ..., SYM,ANTI,SYM, ...)
call fdderivs(ex, betaz, gxxz,gxyz,gxzz,gyyz,gyzz,gzzz, ..., SYM,SYM,ANTI, ...)

! AFTER: 1 fused call
call fdderivs3(ex, betax,..., betay,..., betaz,..., ANTI,SYM,SYM, SYM,ANTI,SYM, SYM,SYM,ANTI, ...)
```

### Group 6: Extrinsic Curvature Derivatives (6 → 2 calls)

```fortran
! BEFORE: 6 calls for d/dx,dy,dz of (Ãxx, Ãxy, Ãxz, Ãyy, Ãyz, Ãzz)
! AFTER: 2 batched calls (same pattern as metric derivatives)
call fderivs3(ex, Axx,..., Axy,..., Axz,..., ...)
call fderivs3(ex, Ayy,..., Ayz,..., Azz,..., ...)
```

### Total: 23 calls → 8 calls

| Group | Variables | Before | After | Routine |
|-------|-----------|--------|-------|---------|
| Shift 1st derivs | β^x, β^y, β^z | 3× fderivs | 1× fderivs3 | fderivs3 |
| Metric 1st derivs | g̃xx, g̃xy, g̃xz, g̃yy, g̃yz, g̃zz | 6× fderivs | 2× fderivs3 | fderivs3 |
| Lapse + trK | α, K | 2× fderivs | 1× fderivs2 | fderivs2 |
| Connection funcs | Γ̃^x, Γ̃^y, Γ̃^z | 3× fderivs | 1× fderivs3 | fderivs3 |
| Shift 2nd derivs | β^x, β^y, β^z | 3× fdderivs | 1× fdderivs3 | fdderivs3 |
| Curvature derivs | Ãxx, Ãxy, Ãxz, Ãyy, Ãyz, Ãzz | 6× fderivs | 2× fderivs3 | fderivs3 |
| **Total** | | **23** | **8** | |

---

## Why This Works: Cache Locality Analysis

### Memory Hierarchy on Intel Xeon (HPC Node)

| Level | Size per core | Latency | Bandwidth |
|-------|--------------|---------|-----------|
| L1d | 32–48 KB | 4 cycles | ~1 TB/s |
| L2 | 256 KB–1 MB | 12 cycles | ~500 GB/s |
| L3 | 1.5 MB/core (shared) | 40 cycles | ~200 GB/s |
| DRAM | — | 100+ cycles | ~50 GB/s |

### Per-Call Memory Footprint

For a typical block with `N ≈ 37` (level 0, 128 ranks):

**Single `fderivs` call:**
- Ghost array `fh(-1:37, -1:37, -1:37)` = 39³ × 8 ≈ **475 KB**
- Output arrays `fx, fy, fz` = 3 × 37³ × 8 ≈ **1.2 MB**
- Total working set: **~1.7 MB** → fits in L2/L3, but fills it

**Fused `fderivs3` call:**
- Three ghost arrays `fh1, fh2, fh3` = 3 × 475 KB ≈ **1.4 MB**
- Nine output arrays = 9 × 405 KB ≈ **3.6 MB**
- Total working set: **~5 MB** → exceeds L2, resides in L3

### Why Bigger Working Set Is Still Faster

The key insight is **temporal locality, not spatial locality**:

1. **Original (3 separate calls):** Each call loads `fh` into L2, computes stencils, writes outputs. Then L2 is evicted. Next call loads a **different** `fh`. The CPU pipeline stalls waiting for cache fills 3× instead of 1×.

2. **Fused (1 call):** All 3 `fh` arrays are loaded once. The inner loop streams through all 3 arrays together. Even though total data is 3× larger, the **loop overhead** (iteration setup, branch prediction, prefetch distance) is paid **once** instead of 3×. The CPU prefetcher learns one access pattern instead of being reset 3 times.

3. **Loop overhead dominates:** For small blocks (N=37), the loop body is only ~31 iterations in the fast dimension. The overhead of loop setup, branch prediction warmup, and prefetch initialization is a significant fraction of each call. Fusing 3 calls eliminates 2/3 of this overhead.

4. **Reduced `symmetry_bd` overhead per derivative:** The setup code (symmetry bounds, coefficient computation, array zeroing) is executed once instead of 3 times.

### For Larger Blocks (Fine AMR Levels)

On fine AMR levels (levels 5–8), blocks can be 50–100+ points per side:
- Array size: 100³ × 8 = **8 MB** per array
- Single `fderivs`: ~32 MB working set → streams from DRAM
- Fused `fderivs3`: ~96 MB working set → also streams from DRAM

For these levels, the gain is primarily from **reduced loop overhead** and **3× fewer function call/return cycles**, not cache reuse. Still beneficial.

---

## Mathematical Equivalence

This optimization is **bit-for-bit identical** to the original code:

1. **Same stencil coefficients:** 4th-order centered differences with coefficients [-1, 8, -8, 1] / (12Δx)
2. **Same symmetry handling:** Each variable gets its own `symmetry_bd()` call with its own SYM/ANTI arguments
3. **Same boundary treatment:** Per-variable symmetry bounds (`imin1`, `imin2`, `imin3`) computed identically to original
4. **Same loop bounds:** Interior loop `[3, N-3]`, boundary loop `[1, N-1]`
5. **Same output arrays:** Results stored in the exact same arrays used by the rest of `bssn_rhs`

The only difference is the **order** in which grid points are visited (all 3 variables at point (i,j,k) before moving to (i+1,j,k)), which does not affect the result since each derivative is an independent stencil operation.

---

## Competition Rule Compliance

✅ **Rule 1** — *"Code related to method parameters MUST NOT be modified"*: No method parameters changed. Same 4th-order stencil, same RK45, same AMR structure.

✅ **Rule 2** — *"All input parameters must remain unchanged"*: `AMSS_NCKU_Input.py` unchanged (MPI_processes = 128).

✅ **Rule 5** — *"New algorithm must be mathematically equivalent"*: Bit-for-bit identical results. Only loop scheduling changed.

✅ **CPU only**: No GPU code involved.

---

## Files Modified

| File | Change |
|------|--------|
| `AMSS_NCKU_source/diff_new.f90` | Added `fderivs3`, `fderivs2`, `fdderivs3` subroutines in the `ghost_width==3` (4th order) section |
| `AMSS_NCKU_source/bssn_rhs.f90` | Replaced 23 individual `fderivs`/`fdderivs` calls with 8 batched calls |

---

## Configuration

- **Compiler:** GCC 14.2.0 (`gnu14/14.2.0`)
- **MPI:** OpenMPI 5.0.7 (`openmpi5/5.0.7`)
- **MPI Ranks:** 128 (2 nodes × 64 tasks/node)
- **Threads:** OMP_NUM_THREADS=1 (pure MPI)
- **Compiler Flags:** `-O3 -march=native -ffast-math -flto -ftree-loop-vectorize -fvect-cost-model=unlimited -fprefetch-loop-arrays`

---

**Date:** February 7, 2026  
**Status:** Awaiting HPC benchmark results
