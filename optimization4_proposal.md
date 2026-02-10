# Optimization 4 — Actual Results

## Summary

This document records the optimization attempts made after commit `d727773` (opt3C baseline: 14.89s/step).
Only **one** of the attempts succeeded. The others were abandoned after testing showed zero improvement or runtime failures.

---

## Current Performance

| Commit | Description | Avg step time | Improvement vs original |
|--------|-------------|---------------|------------------------|
| Original | No optimizations | ~17.56s | — |
| `c290a95` | opt1: compiler flags | ~16.72s | ~4.8% |
| `6f87848` | opt2: loop splitting | ~16.60s | ~5.5% |
| `813a436` | opt2B: fused fderivs3 | ~16.43s | ~6.4% |
| `d727773` | opt3C: kodis + allreduce | 14.89s | ~15.2% |
| **Job 1243** | **opt4: static MPI buffers** | **14.716s** | **~16.2%** |

**Improvement from static buffers alone**: 14.89s → 14.716s = **1.17%**

---

## ✅ Successful: Static MPI Buffers in Parallel.C (Job 1243)

### What was done
Replaced per-call `new`/`delete` of bookkeeping arrays in four MPI communication functions with `static` buffers that are allocated once and reused:

- **`transfer()`**: `MPI_Request*`, `MPI_Status*`, `double**` (send_data, rec_data), `double*` (src)
- **`transfermix()`**: `MPI_Request*`, `MPI_Status*`, `double**` (send_data, rec_data), `double*` (transfer_src, transfer_dst)
- **`Sync(Patch*)`**: `MyList<GridSegment>**` (GC, GO)
- **`Sync(MyList<Patch>*)`**: `MyList<GridSegment>**` (GC, GO)

Pattern used:
```cpp
static T *ptr = nullptr;
static int alloc_size = 0;
if (cpusize > alloc_size) {
    delete[] ptr;
    ptr = new T[cpusize];
    alloc_size = cpusize;
}
```

### Why it works
With ~1022 Sync calls per coarsest timestep (511 Steps × 2 Syncs), each triggering `transfer()` which allocates/frees 5 arrays of 128 elements, this eliminates **~10,000+ malloc/free pairs per coarsest timestep**. The savings come from reduced heap allocator overhead and improved cache behavior.

### Result
- **Job 1243**: avg 14.716s over 23 steps (baseline 14.89s)
- **Improvement**: 1.17%

---

## ❌ Failed: Sync Topology Caching (Jobs 1239, 1240)

### What was attempted
Cache the `GridSegment` lists built during `Sync()` so they don't need to be rebuilt on every call. Since the grid topology doesn't change between Steps, the segment lists (which describe which ghost zones overlap with which blocks on which MPI ranks) should be reusable.

### Implementation
Added a `SyncCacheEntry` struct to `Parallel.h` storing the cached GC/GO segment lists, keyed by the `Patch*` pointer. Used `std::unordered_map<Patch*, SyncCacheEntry>` to look up cached entries.

### What went wrong
**Job 1239**: `MPI_ERR_TRUNCATE` crash during Regrid at timestep 3. Root cause: during Regrid, `Patch*` pointers are freed (`PatL[lev]->destroyList(); PatL[lev] = tmPat;`), and `malloc` reuses the same memory addresses for new Patch objects. The cache then returns stale segment lists for a different grid topology, causing MPI message size mismatches.

**Fix attempt**: Added `invalidate_sync_cache()` calls before Regrid and inside the `recompose_cgh` loop.

**Job 1240**: Still `MPI_ERR_TRUNCATE`. The pointer lifecycle during multi-level Regrid is too complex — new Patch objects are allocated within the Regrid loop itself, and pointer reuse happens before invalidation can run.

### Why it was abandoned
The `Patch*` pointer reuse problem is fundamental to how AMR regridding works. A safe fix would require either (a) a unique generation counter per Patch, or (b) caching at a completely different level (e.g., per-level, keyed by regrid counter). Both are too invasive and risky for a competition setting.

---

## ❌ Failed: Profile-Guided Optimization / PGO (Jobs 1241, 1242)

### What was attempted
3-phase PGO build: (1) compile with `-fprofile-generate`, (2) run a short training simulation, (3) recompile with `-fprofile-use` for optimized production binary.

### What went wrong
**Job 1241** (with `-flto -fprofile-generate`): Instrumented binary was so slow it couldn't complete a single timestep in 9+ hours. The overhead of profiling counters combined with LTO made the binary unusable.

**Job 1242** (without `-flto`, just `-fprofile-generate`): Still stuck for 5+ hours with zero timesteps completed. Root cause: 128 MPI ranks all contending on shared `.gcda` profile data files over the NFS filesystem. GCC's profile instrumentation uses atomic file locking, and with 128 processes on 2 nodes writing to the same NFS-mounted files, the I/O contention completely serialized execution.

### Why it was abandoned
PGO requires either (a) local disk for profile data (not available on this cluster), (b) per-rank profile directories (GCC doesn't support this natively), or (c) running training with far fewer ranks (would produce unrepresentative profiles). None are practical.

---

## ❌ Failed: Opt4 Proposals A, C, D from Original Analysis (Job 1238)

### What was attempted
Three changes tested together:
1. **Remove NaN check** from `compute_rhs_bssn` predictor path
2. **Eliminate redundant `fderivs(chi)`** in GAUGE=2 gauge condition
3. **Skip constraint computation** on predictor calls (`co=0` → `co=1`)

### Result
**Job 1238**: avg 14.922s — **zero improvement** (actually 0.2% worse than baseline 14.89s).

### Why it failed
- The NaN check involves `sum()` over arrays that are already in cache from prior operations — the cost is negligible
- The `fderivs(chi)` call is only 1 out of 84 derivative calls per RHS — 1.2% of derivative work
- The constraint computation in the predictor is cheap compared to the full RHS, and `co=0` is only called once per 4 RK4 substeps
- All three changes combined were within noise

### Status
All three reverted.

---

## Remaining Ideas (Not Yet Attempted)

| Idea | Expected gain | Difficulty | Notes |
|------|--------------|------------|-------|
| Batched symmetry_bd for same-symmetry vars | 2–5% | High | Fuse N symmetry_bd calls into one batch copy. Very invasive to Fortran derivative routines. |
| Memory pool for per-node send/recv arrays | 0.5–1% | Medium | The per-node `double*` arrays in transfer() still do per-call new/delete. Pool allocator could help. |
| Reduce Sync frequency | 3–5% | High | Currently 2 Syncs per Step. If RK4 substeps don't need inter-step Sync, could reduce to 1. Very risky. |
| Loop tiling in compute_rhs_bssn | 1–3% | Very high | Working set ~2.5MB exceeds L2. Tiling could help but extremely invasive. |
