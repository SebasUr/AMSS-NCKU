# Optimization 3 — Kreiss-Oliger Dissipation Loop & MPI\_Allreduce Deferral

**Commit:** `d727773` (`opt3C 1235`)  
**Files modified:** `AMSS_NCKU_source/kodiss.f90`, `AMSS_NCKU_source/bssn_class.C`

---

## Overview

This optimization targets two separate bottlenecks: per-point branch overhead in the Kreiss-Oliger numerical dissipation kernel (`kodis`), and a redundant MPI global synchronization barrier in the RK4 corrector loop.

---

## 1. Kreiss-Oliger Dissipation — Loop Bound Pre-computation (`kodiss.f90`)

### What `kodis` does

The `kodis` subroutine applies 6th-order Kreiss-Oliger dissipation to a single BSSN variable. It uses a 7-point stencil (`i±3`) in each spatial direction. It is called **25 times per RHS evaluation** (once per state variable), and there are **4 RHS evaluations per substep** across **~511 recursive AMR substeps per coarse timestep**. This makes it one of the most frequently invoked kernels.

### The problem

The original code looped over the entire grid (`1..ex`) and checked **at every point** whether the 7-point stencil fit within bounds:

```fortran
do k=1,ex(3)
do j=1,ex(2)
do i=1,ex(1)
  if(i-3 >= imin .and. i+3 <= imax .and. &
     j-3 >= jmin .and. j+3 <= jmax .and. &
     k-3 >= kmin .and. k+3 <= kmax) then
     ! ... compute dissipation ...
  endif
enddo; enddo; enddo
```

This means every grid point evaluates 6 comparisons and a compound boolean before doing any useful work. For interior points (the vast majority), the check always passes. This creates:

- **Branch misprediction overhead** on every iteration.
- **Prevents vectorization** — the compiler cannot guarantee the loop body executes on contiguous iterations, breaking SIMD auto-vectorization.

### The fix

Pre-compute the valid loop bounds *once* before entering the loop:

```fortran
ilo = max(1, imin + 3)
ihi = min(ex(1), imax - 3)
jlo = max(1, jmin + 3)
jhi = min(ex(2), jmax - 3)
klo = max(1, kmin + 3)
khi = min(ex(3), kmax - 3)

if(ilo > ihi .or. jlo > jhi .or. klo > khi) return   ! early exit

do k=klo,khi
do j=jlo,jhi
do i=ilo,ihi
  ! ... compute dissipation (no branch) ...
enddo; enddo; enddo
```

### Why it works

1. **Eliminates per-point branching.** The 6-comparison `if` block is replaced by fixed loop bounds computed once. The inner loop body executes unconditionally for every iteration — zero branch overhead.

2. **Enables auto-vectorization.** With known contiguous bounds and no conditional, the compiler can vectorize the inner `i`-loop with AVX2/AVX-512 SIMD instructions. The stencil arithmetic (additions, multiplications by constants) is perfectly suited for SIMD.

3. **Early exit for empty blocks.** If a block has no valid interior points (e.g., very small ghost-only patches), the subroutine returns immediately without even calling `symmetry_bd`.

### Additional micro-optimizations

- **Division → multiplication:** `eps/cof` and `1/dX`, `1/dY`, `1/dZ` are pre-computed as `ecof`, `rdX`, `rdY`, `rdZ` outside the loop. Inside the loop, `/dX` becomes `*rdX`. Division is 4–10× slower than multiplication on modern CPUs, and this eliminates 3 divisions per grid point.

- **Dead code removal:** The `#if 0 / #else` preprocessor block (separate x/y/z passes vs. combined) was cleaned up, keeping only the combined version that performs one write per point instead of three.

---

## 2. MPI\_Allreduce Deferral in RK4 Corrector (`bssn_class.C`)

### What it did

In the `Step()` function, each RK4 substep performed:
1. Compute RHS
2. Apply RK4 update + Sommerfeld boundary conditions
3. **`MPI_Allreduce` to check for NaN across all 128 ranks**
4. Synchronize ghost zones

This NaN check existed in both the **predictor** (substep 0) and all three **corrector** iterations (substeps 1–3), totaling **4 global reductions per Step() call**.

### The problem

`MPI_Allreduce` is a **global barrier** — every rank must reach the call and exchange data before any can proceed. With 128 ranks across 2 nodes, this introduces:

- **Synchronization latency** (~5–20 µs per call from inter-node fabric).
- **Load imbalance amplification** — the slowest rank determines when all ranks proceed.

With ~511 `Step()` calls per coarse timestep and 3 corrector iterations each, that's **~1,533 unnecessary `MPI_Allreduce` calls per coarse timestep**.

### The fix

Remove the `MPI_Allreduce` NaN check from the corrector iterations, keeping only the predictor check:

```cpp
// NaN check deferred to predictor step of next timestep to avoid
// MPI_Allreduce barrier overhead in corrector iterations
```

### Why it's safe

- If a NaN appears during a corrector substep, it will **propagate** and be caught by the predictor NaN check at the start of the **next** timestep's `Step()` call.  
- NaN arithmetic is "sticky" — any operation involving NaN produces NaN, so it cannot silently disappear.  
- The predictor check (which is kept) still calls `MPI_Allreduce`, detects the NaN, dumps diagnostic data, and aborts via `MPI_Abort`.  
- The only difference is detection is delayed by at most 3 substeps — the simulation still aborts cleanly with full diagnostics.

### Why it helps

Removing 3 out of 4 `MPI_Allreduce` calls per `Step()` eliminates **75% of the NaN-check synchronization barriers**. Each barrier costs synchronization latency plus load-imbalance wait time. Over ~511 calls per coarse timestep, this adds up to measurable wall-clock savings, especially at scale where inter-node latency dominates.

---

## Impact

Both changes target high-frequency operations (called thousands of times per coarse timestep) with low per-call overhead that accumulates significantly:

| Change | Calls removed/optimized per coarse timestep | Mechanism |
|--------|----------------------------------------------|-----------|
| kodis loop bounds | ~51,100 calls × branch elimination | Vectorization + no branch |
| MPI\_Allreduce deferral | ~1,533 global reductions removed | Reduced synchronization |

Combined with the prior optimizations (compiler flags, loop splitting, fused derivatives), these changes contribute to the cumulative performance improvement.
